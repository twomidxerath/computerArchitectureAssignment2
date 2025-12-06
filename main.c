#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <limits.h>

#define NUM_CACHE 5
#define NUM_BLOCK 5
#define NUM_ASSOC 4
#define NUM_ROWS 8
#define NUM_COLS 25
#define MAX_SETS 2048
#define MAX_WAYS 8

static const int CACHE_SIZES[NUM_CACHE] = {1024, 2048, 4096, 8192, 16384};
static const int BLOCK_SIZES[NUM_BLOCK] = {8, 16, 32, 64, 128};
static const int ASSOC_LIST[NUM_ASSOC] = {1, 2, 4, 8};


/* LRU */
struct Block_LRU {
    unsigned long tag[MAX_WAYS];      
    char valid[MAX_WAYS];             
    char write_back[MAX_WAYS];   
    unsigned long lru_time[MAX_WAYS]; 
};
static struct Block_LRU icah_lru[MAX_SETS];
static struct Block_LRU dcah_lru[MAX_SETS];

/* FIFO */
struct Block_FIFO {
    unsigned long tag[MAX_WAYS]; /* [수정됨] char -> unsigned long */
    char valid[MAX_WAYS];
    char write_back[MAX_WAYS];
};
static struct Block_FIFO icah_fifo[MAX_SETS];
static struct Block_FIFO dcah_fifo[MAX_SETS];
static int icah_fifo_ptr[MAX_SETS];
static int dcah_fifo_ptr[MAX_SETS];

static int get_log(int val);
static void calculate_address_fields(unsigned long addr, int block_size, int assoc, int cache_size, unsigned long *tag_out, unsigned long *index_out);

static void simulate_lru  (int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]);

static void simulate_fifo (int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]);

static void simulate_new  (int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]);

static void print_results(const char *label,
                        const double miss[NUM_ROWS][NUM_COLS],
                        const int writes[NUM_ROWS][NUM_COLS]);

static void print_best_results(
    const double lru_miss[NUM_ROWS][NUM_COLS], const int lru_writes[NUM_ROWS][NUM_COLS],
    const int lru_i_tot[NUM_ROWS][NUM_COLS], const int lru_d_tot[NUM_ROWS][NUM_COLS],
    const double fifo_miss[NUM_ROWS][NUM_COLS], const int fifo_writes[NUM_ROWS][NUM_COLS],
    const int fifo_i_tot[NUM_ROWS][NUM_COLS], const int fifo_d_tot[NUM_ROWS][NUM_COLS],
    int i_hit, int i_miss, int d_hit, int d_miss);

static void read_trace(const char *path,
                        int **ptype, unsigned long **paddr, int *plen);

static void print_two_cache_state(const char *policy, 
                        int is_icache, int index, int assoc);

static void print_new_cache_state(int is_icache, int index, int assoc);

static void usage(const char *prog) {
    fprintf(stderr,
            "Usage: %s <policy> <trace_file> [cycle_params]\n"
            "  <policy>        FIFO, LRU, NEW or BEST (case-insensitive)\n"
            "  <trace_file>    input trace in .din format\n"
            "  [cycle_params]  Required only for BEST policy:\n"
            "                    <i_hit> <i_miss> <d_hit> <d_miss>\n"
            "  Example (FIFO):  %s FIFO trace.din\n"
            "  Example (LRU):   %s LRU trace.din\n"
            "  Example (NEW):   %s NEW trace.din\n"
            "  Example (BEST):  %s BEST trace.din 1 100 1 50\n",
            prog, prog, prog, prog, prog);
    exit(1);
}

static void read_trace(const char *path, int **ptype, unsigned long **paddr, int *plen) {
    FILE *fp;
    fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Error opening file %s\n", path); exit(1); }
    
    int count = 0;
    int tmp_time, tmp_label;
    unsigned long tmp_addr;
    
    // 첫 번째 패스: 라인 수 세기
    while(fscanf(fp, "%d %d %lx", &tmp_time, &tmp_label, &tmp_addr) != EOF) {
        count++;
    }
    
    rewind(fp);
    *ptype = (int *)malloc(sizeof(int) * count);
    *paddr = (unsigned long *)malloc(sizeof(unsigned long) * count);

    // 두 번째 패스: 데이터 저장
    for(int i = 0; i < count; i++) {
        fscanf(fp, "%d %d %lx", &tmp_time, &(*ptype)[i], &(*paddr)[i]);
    }
    
    *plen = count;
    fclose(fp);
}


static int get_log(int val){
    int count = 0;
    while(val > 1){
        val >>= 1;
        count++;
    }
    return count;
}
static void calculate_address_fields(unsigned long addr, int block_size, int assoc, int cache_size, unsigned long *tag_out, unsigned long *index_out){
    int set_num = cache_size / (block_size * assoc);
    int offset_bits = get_log(block_size);
    int index_bits = get_log(set_num);
    *index_out = (addr >> offset_bits) & (set_num - 1);
    *tag_out = addr >> (offset_bits + index_bits);
}

static void simulate_lru(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {

    int i, k, l;
    int assoc_idx, block_idx, cache_idx;
    unsigned long current_time = 0;

    for (assoc_idx = 0; assoc_idx < NUM_ASSOC; assoc_idx++) {
        int assoc = ASSOC_LIST[assoc_idx];

        for (block_idx = 0; block_idx < NUM_BLOCK; block_idx++) {
            int blk_size = BLOCK_SIZES[block_idx];

            for (cache_idx = 0; cache_idx < NUM_CACHE; cache_idx++) {
                int cache_size = CACHE_SIZES[cache_idx];

                memset(icah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                memset(dcah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                
                int i_miss_count = 0, d_miss_count = 0;
                int i_acc_count = 0, d_acc_count = 0;
                int d_write_mem_count = 0;

                for (k = 0; k < length; k++) {
                    current_time++;
                    int cur_type = type[k];
                    unsigned long cur_addr = addr[k];
                    unsigned long tag_val, set_index;
                    
                    calculate_address_fields(cur_addr, blk_size, assoc, cache_size, &tag_val, &set_index);

                    struct Block_LRU *target_cache = (cur_type == 2) ? icah_lru : dcah_lru;
                    if (cur_type == 2) i_acc_count++; else d_acc_count++;

                    int hit = 0, hit_way = -1;
                    for (l = 0; l < assoc; l++) {
                        if (target_cache[set_index].valid[l] && target_cache[set_index].tag[l] == tag_val) {
                            hit = 1; hit_way = l;
                            if (cur_type == 1) target_cache[set_index].write_back[l] = 1;
                            break; 
                        }
                    }

                    if (hit) {
                        target_cache[set_index].lru_time[hit_way] = current_time;
                    } else {
                        if (cur_type == 2) i_miss_count++; else d_miss_count++;

                        int victim_way = -1;
                        unsigned long oldest_time = ULONG_MAX; // limits.h 필요
                        int found_empty = 0;

                        for (l = 0; l < assoc; l++) {
                            if (target_cache[set_index].valid[l] == 0) {
                                victim_way = l; found_empty = 1; break;
                            }
                        }

                        if (!found_empty) {
                            for (l = 0; l < assoc; l++) {
                                if (target_cache[set_index].lru_time[l] < oldest_time) {
                                    oldest_time = target_cache[set_index].lru_time[l];
                                    victim_way = l;
                                }
                            }
                        }

                        if (target_cache[set_index].valid[victim_way] && target_cache[set_index].write_back[victim_way]) {
                            d_write_mem_count++;
                        }

                        target_cache[set_index].valid[victim_way] = 1;
                        target_cache[set_index].tag[victim_way] = tag_val;
                        target_cache[set_index].lru_time[victim_way] = current_time;
                        target_cache[set_index].write_back[victim_way] = (cur_type == 1) ? 1 : 0;
                    }
                }

                int row_i = assoc_idx;
                int row_d = assoc_idx + NUM_ASSOC;
                int col = (block_idx * NUM_CACHE) + cache_idx;

                miss[row_i][col] = (i_acc_count > 0) ? (double)i_miss_count / i_acc_count : 0.0;
                miss[row_d][col] = (d_acc_count > 0) ? (double)d_miss_count / d_acc_count : 0.0;
                writes[row_d][col] = d_write_mem_count;
                writes[row_i][col] = 0;
                i_totals[row_i][col] = i_acc_count;
                d_totals[row_d][col] = d_acc_count;
            }
        }
    }
}

static void simulate_fifo(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {

    int i, k, l;
    int assoc_idx, block_idx, cache_idx;

    for (assoc_idx = 0; assoc_idx < NUM_ASSOC; assoc_idx++) {
        int assoc = ASSOC_LIST[assoc_idx];
        for (block_idx = 0; block_idx < NUM_BLOCK; block_idx++) {
            int blk_size = BLOCK_SIZES[block_idx];
            for (cache_idx = 0; cache_idx < NUM_CACHE; cache_idx++) {
                int cache_size = CACHE_SIZES[cache_idx];

                memset(icah_fifo, 0, sizeof(struct Block_FIFO) * MAX_SETS);
                memset(dcah_fifo, 0, sizeof(struct Block_FIFO) * MAX_SETS);
                memset(icah_fifo_ptr, 0, sizeof(int) * MAX_SETS);
                memset(dcah_fifo_ptr, 0, sizeof(int) * MAX_SETS);

                int i_miss_count = 0, d_miss_count = 0;
                int i_acc_count = 0, d_acc_count = 0;
                int d_write_mem_count = 0;

                for (k = 0; k < length; k++) {
                    int cur_type = type[k];
                    unsigned long cur_addr = addr[k];
                    unsigned long tag_val, set_index;
                    
                    calculate_address_fields(cur_addr, blk_size, assoc, cache_size, &tag_val, &set_index);

                    struct Block_FIFO *target_cache;
                    int *target_ptr;
                    if (cur_type == 2) { target_cache = icah_fifo; target_ptr = icah_fifo_ptr; i_acc_count++; }
                    else { target_cache = dcah_fifo; target_ptr = dcah_fifo_ptr; d_acc_count++; }

                    int hit = 0;
                    for (l = 0; l < assoc; l++) {
                        if (target_cache[set_index].valid[l] && target_cache[set_index].tag[l] == tag_val) {
                            hit = 1;
                            if (cur_type == 1) target_cache[set_index].write_back[l] = 1;
                            break; 
                        }
                    }

                    if (!hit) {
                        if (cur_type == 2) i_miss_count++; else d_miss_count++;
                        
                        int victim_way = target_ptr[set_index];
                        if (target_cache[set_index].valid[victim_way] && target_cache[set_index].write_back[victim_way]) {
                            d_write_mem_count++;
                        }

                        target_cache[set_index].valid[victim_way] = 1;
                        target_cache[set_index].tag[victim_way] = tag_val;
                        target_cache[set_index].write_back[victim_way] = (cur_type == 1) ? 1 : 0;
                        
                        target_ptr[set_index] = (victim_way + 1) % assoc;
                    }
                }

                int row_i = assoc_idx;
                int row_d = assoc_idx + NUM_ASSOC;
                int col = (block_idx * NUM_CACHE) + cache_idx;

                miss[row_i][col] = (i_acc_count > 0) ? (double)i_miss_count / i_acc_count : 0.0;
                miss[row_d][col] = (d_acc_count > 0) ? (double)d_miss_count / d_acc_count : 0.0;
                writes[row_d][col] = d_write_mem_count;
                writes[row_i][col] = 0;
                i_totals[row_i][col] = i_acc_count;
                d_totals[row_d][col] = d_acc_count;
            }
        }
    }
}

static void simulate_new(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {
    // MRU Implementation
    int i, k, l;
    int assoc_idx, block_idx, cache_idx;
    unsigned long current_time = 0;

    for (assoc_idx = 0; assoc_idx < NUM_ASSOC; assoc_idx++) {
        int assoc = ASSOC_LIST[assoc_idx];
        for (block_idx = 0; block_idx < NUM_BLOCK; block_idx++) {
            int blk_size = BLOCK_SIZES[block_idx];
            for (cache_idx = 0; cache_idx < NUM_CACHE; cache_idx++) {
                int cache_size = CACHE_SIZES[cache_idx];

                memset(icah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                memset(dcah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                
                int i_miss_count = 0, d_miss_count = 0;
                int i_acc_count = 0, d_acc_count = 0;
                int d_write_mem_count = 0;

                for (k = 0; k < length; k++) {
                    current_time++;
                    int cur_type = type[k];
                    unsigned long cur_addr = addr[k];
                    unsigned long tag_val, set_index;
                    
                    calculate_address_fields(cur_addr, blk_size, assoc, cache_size, &tag_val, &set_index);

                    struct Block_LRU *target_cache = (cur_type == 2) ? icah_lru : dcah_lru;
                    if (cur_type == 2) i_acc_count++; else d_acc_count++;

                    int hit = 0, hit_way = -1;
                    for (l = 0; l < assoc; l++) {
                        if (target_cache[set_index].valid[l] && target_cache[set_index].tag[l] == tag_val) {
                            hit = 1; hit_way = l;
                            if (cur_type == 1) target_cache[set_index].write_back[l] = 1;
                            break; 
                        }
                    }

                    if (hit) {
                        target_cache[set_index].lru_time[hit_way] = current_time;
                    } else {
                        if (cur_type == 2) i_miss_count++; else d_miss_count++;

                        int victim_way = -1;
                        unsigned long most_recent_time = 0; // 0으로 초기화 (MRU)
                        int found_empty = 0;

                        for (l = 0; l < assoc; l++) {
                            if (target_cache[set_index].valid[l] == 0) {
                                victim_way = l; found_empty = 1; break;
                            }
                        }

                        if (!found_empty) {
                            for (l = 0; l < assoc; l++) {
                                // MRU: 가장 시간값이 큰(최근) 블록 찾기
                                if (target_cache[set_index].lru_time[l] > most_recent_time) {
                                    most_recent_time = target_cache[set_index].lru_time[l];
                                    victim_way = l;
                                }
                            }
                        }

                        if (target_cache[set_index].valid[victim_way] && target_cache[set_index].write_back[victim_way]) {
                            d_write_mem_count++;
                        }

                        target_cache[set_index].valid[victim_way] = 1;
                        target_cache[set_index].tag[victim_way] = tag_val;
                        target_cache[set_index].lru_time[victim_way] = current_time;
                        target_cache[set_index].write_back[victim_way] = (cur_type == 1) ? 1 : 0;
                    }
                    
                    // 디버깅용 (보고서 제출 시 활용)
                    if (k < 10) print_new_cache_state((cur_type==2), set_index, assoc);
                }

                int row_i = assoc_idx;
                int row_d = assoc_idx + NUM_ASSOC;
                int col = (block_idx * NUM_CACHE) + cache_idx;

                miss[row_i][col] = (i_acc_count > 0) ? (double)i_miss_count / i_acc_count : 0.0;
                miss[row_d][col] = (d_acc_count > 0) ? (double)d_miss_count / d_acc_count : 0.0;
                writes[row_d][col] = d_write_mem_count;
                writes[row_i][col] = 0;
                i_totals[row_i][col] = i_acc_count;
                d_totals[row_d][col] = d_acc_count;
            }
        }
    }
}

static void print_results(const char *label, const double miss[NUM_ROWS][NUM_COLS], const int writes[NUM_ROWS][NUM_COLS]) {
    int i, j, k;
    printf("\nMissRate\n");
    for (i = 0; i < NUM_ROWS; i++) {
        if (i == 0) {
            printf("           ");
            for (k = 0; k < NUM_BLOCK; k++) printf("%s/%-4d                                 ", label, BLOCK_SIZES[k]);
            printf("\nI cache   ");
            for (k = 0; k < NUM_BLOCK; k++) for (j = 0; j < NUM_CACHE; j++) printf("%-7d", CACHE_SIZES[j]);
            printf("\n");
        }
        if (i == NUM_ASSOC) {
            printf("\n           ");
            for (k = 0; k < NUM_BLOCK; k++) printf("%s/%-4d                                 ", label, BLOCK_SIZES[k]);
            printf("\nD cache   ");
            for (k = 0; k < NUM_BLOCK; k++) for (j = 0; j < NUM_CACHE; j++) printf("%-7d", CACHE_SIZES[j]);
            printf("\n");
        }
        if (i % NUM_ASSOC == 0) printf("Direct | ");
        else if (i % NUM_ASSOC == 1) printf("2  Way | ");
        else if (i % NUM_ASSOC == 2) printf("4  Way | ");
        else printf("8  Way | ");
        
        for (j = 0; j < NUM_COLS; j++) printf("%.4lf ", miss[i][j]);
        printf("\n");
    }

    printf("\nWrite Count\n");
    for (i = 0; i < NUM_ROWS; i++) {
        if (i == 0) {
            printf("           ");
            for (k = 0; k < NUM_BLOCK; k++) printf("%s/%-4d                          ", label, BLOCK_SIZES[k]);
            printf("\nI cache   ");
            for (k = 0; k < NUM_BLOCK; k++) for (j = 0; j < NUM_CACHE; j++) printf("%6d", CACHE_SIZES[j]);
            printf("\n");
        }
        if (i == NUM_ASSOC) {
            printf("\n           ");
            for (k = 0; k < NUM_BLOCK; k++) printf("%s/%-4d                          ", label, BLOCK_SIZES[k]);
            printf("\nD cache   ");
            for (k = 0; k < NUM_BLOCK; k++) for (j = 0; j < NUM_CACHE; j++) printf("%6d", CACHE_SIZES[j]);
            printf("\n");
        }
        if (i % NUM_ASSOC == 0) printf("Direct | ");
        else if (i % NUM_ASSOC == 1) printf("2  Way | ");
        else if (i % NUM_ASSOC == 2) printf("4  Way | ");
        else printf("8  Way | ");

        for (j = 0; j < NUM_COLS; j++) printf("%5d ", writes[i][j]);
        printf("\n");
    }
}

static void print_best_results(
    const double lru_miss[NUM_ROWS][NUM_COLS], const int lru_writes[NUM_ROWS][NUM_COLS],
    const int lru_i_tot[NUM_ROWS][NUM_COLS], const int lru_d_tot[NUM_ROWS][NUM_COLS],
    const double fifo_miss[NUM_ROWS][NUM_COLS], const int fifo_writes[NUM_ROWS][NUM_COLS],
    const int fifo_i_tot[NUM_ROWS][NUM_COLS], const int fifo_d_tot[NUM_ROWS][NUM_COLS],
    int i_hit, int i_miss, int d_hit, int d_miss) {
    
    int cl, bl, as; 

    // 1. 캐시 크기별로 최적의 설정을 찾습니다.
    for (cl = 0; cl < NUM_CACHE; cl++) {

        // [핵심 수정] 변수를 2개로 분리했습니다!
        double min_i_cycles = DBL_MAX;
        double min_d_cycles = DBL_MAX;

        // 최적 설정 저장 변수들도 분리
        const char *best_i_policy = "N/A";
        int best_i_block = 0;
        int best_i_assoc = 0;
        double best_i_missrate = 0.0;
        double best_i_total_cycles = 0.0;

        const char *best_d_policy = "N/A";
        int best_d_block = 0;
        int best_d_assoc = 0;
        double best_d_missrate = 0.0;
        int best_d_writes = 0;
        double best_d_total_cycles = 0.0;

        // 2. 모든 블록 크기와 연관도 조합을 순회
        for (bl = 0; bl < NUM_BLOCK; bl++) {
            for (as = 0; as < NUM_ASSOC; as++) {
                
                int row_i = as;
                int row_d = as + NUM_ASSOC;
                int col = (bl * NUM_CACHE) + cl;
                
                // --- (A) LRU 성능 계산 ---
                double lru_i_mr = lru_miss[row_i][col];
                int lru_i_total = lru_i_tot[row_i][col];
                // Hit Count = Total * (1 - MissRate)
                double lru_i_cycles = (lru_i_total * (1.0 - lru_i_mr) * i_hit) + 
                                      (lru_i_total * lru_i_mr * i_miss);

                // [수정] I-Cache만 따로 비교 (LRU)
                if (lru_i_cycles < min_i_cycles) {
                    min_i_cycles = lru_i_cycles;
                    best_i_policy = "LRU";
                    best_i_block = BLOCK_SIZES[bl];
                    best_i_assoc = ASSOC_LIST[as];
                    best_i_missrate = lru_i_mr;
                    best_i_total_cycles = lru_i_cycles;
                }

                double lru_d_mr = lru_miss[row_d][col];
                int lru_d_total = lru_d_tot[row_d][col];
                double lru_d_cycles = (lru_d_total * (1.0 - lru_d_mr) * d_hit) + 
                                      (lru_d_total * lru_d_mr * d_miss);

                // [수정] D-Cache만 따로 비교 (LRU)
                if (lru_d_cycles < min_d_cycles) {
                    min_d_cycles = lru_d_cycles;
                    best_d_policy = "LRU";
                    best_d_block = BLOCK_SIZES[bl];
                    best_d_assoc = ASSOC_LIST[as];
                    best_d_missrate = lru_d_mr;
                    best_d_writes = lru_writes[row_d][col];
                    best_d_total_cycles = lru_d_cycles;
                }

                // --- (B) FIFO 성능 계산 ---
                double fifo_i_mr = fifo_miss[row_i][col];
                int fifo_i_total = fifo_i_tot[row_i][col];
                double fifo_i_cycles = (fifo_i_total * (1.0 - fifo_i_mr) * i_hit) + 
                                       (fifo_i_total * fifo_i_mr * i_miss);

                // [수정] I-Cache만 따로 비교 (FIFO)
                if (fifo_i_cycles < min_i_cycles) {
                    min_i_cycles = fifo_i_cycles;
                    best_i_policy = "FIFO";
                    best_i_block = BLOCK_SIZES[bl];
                    best_i_assoc = ASSOC_LIST[as];
                    best_i_missrate = fifo_i_mr;
                    best_i_total_cycles = fifo_i_cycles;
                }

                double fifo_d_mr = fifo_miss[row_d][col];
                int fifo_d_total = fifo_d_tot[row_d][col];
                double fifo_d_cycles = (fifo_d_total * (1.0 - fifo_d_mr) * d_hit) + 
                                       (fifo_d_total * fifo_d_mr * d_miss);

                // [수정] D-Cache만 따로 비교 (FIFO)
                if (fifo_d_cycles < min_d_cycles) {
                    min_d_cycles = fifo_d_cycles;
                    best_d_policy = "FIFO";
                    best_d_block = BLOCK_SIZES[bl];
                    best_d_assoc = ASSOC_LIST[as];
                    best_d_missrate = fifo_d_mr;
                    best_d_writes = fifo_writes[row_d][col];
                    best_d_total_cycles = fifo_d_cycles;
                }
            }
        }

        // 3. 결과 출력
        printf("--- Cache Size: %d bytes ---\n", CACHE_SIZES[cl]);

        if (min_i_cycles == DBL_MAX)
            printf("  I-Cache: No instruction accesses.\n");
        else
            printf("  Best I-Cache: Policy=%-4s | Block=%-4d | Assoc=%-2d | MissRate=%.4f | Total Cycles=%.0f\n",
                   best_i_policy, best_i_block, best_i_assoc, best_i_missrate, best_i_total_cycles);

        if (min_d_cycles == DBL_MAX)
            printf("  D-Cache: No data accesses.\n");
        else
            printf("  Best D-Cache: Policy=%-4s | Block=%-4d | Assoc=%-2d | MissRate=%.4f | Writes=%-5d | Total Cycles=%.0f\n",
                   best_d_policy, best_d_block, best_d_assoc, best_d_missrate, best_d_writes, best_d_total_cycles);

        printf("\n");
    }
}

static void print_two_cache_state(const char *policy, int is_icache, int index, int assoc) {
    printf("\n[Cache State Dump] Policy=%s | %s | index=%d | assoc=%d\n", policy, is_icache ? "I-Cache" : "D-Cache", index, assoc);
    if (strcmp(policy, "LRU") == 0 || strcmp(policy, "lru") == 0) {
        struct Block_LRU *cache = is_icache ? icah_lru : dcah_lru;
        for (int i = 0; i < assoc; i++) {
            printf("  Way %-2d | valid=%d  tag=%lu  write_back=%d\n", i, cache[index].valid[i], cache[index].tag[i], cache[index].write_back[i]);
        }
    } else if (strcmp(policy, "FIFO") == 0 || strcmp(policy, "fifo") == 0) {
        struct Block_FIFO *cache = is_icache ? icah_fifo : dcah_fifo;
        int *ptr = is_icache ? icah_fifo_ptr : dcah_fifo_ptr;
        printf("  FIFO pointer = %d\n", ptr[index]);
        for (int i = 0; i < assoc; i++) {
            printf("  Way %-2d | valid=%d  tag=%lu  write_back=%d\n", i, cache[index].valid[i], cache[index].tag[i], cache[index].write_back[i]);
        }
    }
    printf("\n");
}

static void print_new_cache_state(int is_icache, int index, int assoc) {
    printf("\n[NEW Policy State Dump] %s | index=%d | assoc=%d\n", is_icache ? "I-Cache" : "D-Cache", index, assoc);
    struct Block_LRU *cache = is_icache ? icah_lru : dcah_lru;
    for (int i = 0; i < assoc; i++) {
        printf("  Way %-2d | valid=%d  tag=%lx  write_back=%d  lru_time=%lu\n", i, cache[index].valid[i], cache[index].tag[i], cache[index].write_back[i], cache[index].lru_time[i]); 
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    if (argc < 3 || (argc > 3 && argc < 7) || argc > 7) usage(argv[0]);

    int policy = 0;  
    int i_hit_c = 0, i_miss_c = 0, d_hit_c = 0, d_miss_c = 0;
    char *trace_file = NULL;

    if (!strcasecmp(argv[1], "FIFO")) { if (argc != 3) usage(argv[0]); policy = 1; trace_file = argv[2]; }
    else if (!strcasecmp(argv[1], "LRU")) { if (argc != 3) usage(argv[0]); policy = 0; trace_file = argv[2]; }
    else if (!strcasecmp(argv[1], "NEW")) { if (argc != 3) usage(argv[0]); policy = 3; trace_file = argv[2]; }
    else if (!strcasecmp(argv[1], "BEST")) { if (argc != 7) usage(argv[0]); policy = 2; trace_file = argv[2]; i_hit_c = atoi(argv[3]); i_miss_c = atoi(argv[4]); d_hit_c = atoi(argv[5]); d_miss_c = atoi(argv[6]); }
    else { usage(argv[0]); }

    int *type = NULL;
    unsigned long *addr = NULL;
    int length = 0;

    printf("Reading trace file: %s\n", trace_file); // [수정 2] 정답 파일(sample_output) 형식에 맞춰 전체 경로 출력 유지

    read_trace(trace_file, &type, &addr, &length);
    printf("Trace contains %d memory accesses.\n", length);

    if (policy == 0) {  
        printf("Simulating LRU policy...\n");
        double miss[NUM_ROWS][NUM_COLS] = {{0}};
        int writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int d_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        simulate_lru(type, addr, length, miss, writes, i_tot, d_tot);
        print_results("LRU", miss, writes);
    } else if (policy == 1) { 
        printf("Simulating FIFO policy...\n");
        double miss[NUM_ROWS][NUM_COLS] = {{0}};
        int writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int d_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        simulate_fifo(type, addr, length, miss, writes, i_tot, d_tot);
        print_results("FIFO", miss, writes);
    } else if (policy == 3) { 
        printf("Simulating NEW policy...\n");
        double miss[NUM_ROWS][NUM_COLS] = {{0}};
        int writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int d_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        simulate_new(type, addr, length, miss, writes, i_tot, d_tot);
        print_results("NEW", miss, writes);
    } else { 
        printf("Simulating LRU policy for BEST...\n");
        double lru_miss[NUM_ROWS][NUM_COLS] = {{0}};
        int lru_writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int lru_i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int lru_d_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        simulate_lru(type, addr, length, lru_miss, lru_writes, lru_i_tot, lru_d_tot);

        double fifo_miss[NUM_ROWS][NUM_COLS] = {{0}};
        int fifo_writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int fifo_i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int fifo_d_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        printf("Simulating FIFO policy for BEST...\n");
        simulate_fifo(type, addr, length, fifo_miss, fifo_writes, fifo_i_tot, fifo_d_tot);

        printf("\n--- BEST Configuration Analysis ---\n");
        printf("Cycle Parameters: I(Hit/Miss) = %d/%d, D(Hit/Miss) = %d/%d\n\n", i_hit_c, i_miss_c, d_hit_c, d_miss_c);
        print_best_results(lru_miss, lru_writes, lru_i_tot, lru_d_tot, fifo_miss, fifo_writes, fifo_i_tot, fifo_d_tot, i_hit_c, i_miss_c, d_hit_c, d_miss_c);
    }

    free(type);
    free(addr);
    return 0;
}