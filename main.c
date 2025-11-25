#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

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
    char tag[8];
    char valid[8];
    char write_back[8];
};
static struct Block_FIFO icah_fifo[MAX_SETS];
static struct Block_FIFO dcah_fifo[MAX_SETS];
static int icah_fifo_ptr[MAX_SETS];
static int dcah_fifo_ptr[MAX_SETS];

static int log(int val);
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
    int count = 0;
    int tmp_time, tmp_label;
    unsigned long tmp_addr;
    for(count = 0; fscanf(fp, "%d %d %lx", &tmp_time, &tmp_label, &tmp_addr) != EOF; count++);
    rewind(fp);
    *ptype = (int *)malloc(sizeof(int) * count);
    *paddr = (unsigned long *)malloc(sizeof(unsigned long) * count);

    for(int i = 0; fscanf(fp, "%d %d %lx", &tmp_time, &(*ptype)[i], &(*paddr)[i]) != EOF; i++);
    *plen = count;
    fclose(fp);

}


static int log(int val){
    int count = 0;
    while(val > 0){
        val >>= 1;
        count++;
    }
    return count;
}
static void calculate_address_fields(unsigned long addr, int block_size, int assoc, int cache_size, unsigned long *tag_out, unsigned long *index_out){
    int set_num = cache_size / (block_size * assoc);
    int offset_bits = log(block_size);
    int index_bits = log(set_num);
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
    
    /* 시간 추적을 위한 변수 (시뮬레이션 전체 시간) */
    unsigned long current_time = 0;

    /* 3중 루프: 모든 구성 시뮬레이션 */
    for (assoc_idx = 0; assoc_idx < NUM_ASSOC; assoc_idx++) {
        int assoc = ASSOC_LIST[assoc_idx];

        for (block_idx = 0; block_idx < NUM_BLOCK; block_idx++) {
            int blk_size = BLOCK_SIZES[block_idx];

            for (cache_idx = 0; cache_idx < NUM_CACHE; cache_idx++) {
                int cache_size = CACHE_SIZES[cache_idx];

                /* ------------------------------------------------ */
                /* 1. 초기화 (Initialization) */
                /* ------------------------------------------------ */
                
                // 캐시 메모리 및 포인터 초기화
                memset(icah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                memset(dcah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                
                // 통계 변수 초기화
                int i_miss_count = 0;
                int d_miss_count = 0;
                int i_acc_count = 0;
                int d_acc_count = 0;
                int d_write_mem_count = 0;

                /* ------------------------------------------------ */
                /* 2. 트레이스 시뮬레이션 (Trace Loop) */
                /* ------------------------------------------------ */
                for (k = 0; k < length; k++) {
                    current_time++; // ⭐ 시간 흐름 (매 접근마다 증가)

                    int cur_type = type[k];
                    unsigned long cur_addr = addr[k];

                    // 주소 디코딩 (헬퍼 함수)
                    unsigned long tag_val;
                    unsigned long set_index;
                    calculate_address_fields(cur_addr, blk_size, assoc, cache_size, 
                                             &tag_val, &set_index);

                    // 타겟 캐시 설정
                    struct Block_LRU *target_cache;
                    int is_icache = (cur_type == 2);

                    if (is_icache) {
                        target_cache = icah_lru;
                        i_acc_count++;
                    } else {
                        target_cache = dcah_lru;
                        d_acc_count++;
                    }

                    // Hit 검사
                    int hit = 0;
                    int hit_way = -1;

                    for (l = 0; l < assoc; l++) {
                        if (target_cache[set_index].valid[l] && target_cache[set_index].tag[l] == tag_val) {
                            hit = 1;
                            hit_way = l;
                            
                            // Hit & Write -> Dirty Bit 설정
                            if (!is_icache && cur_type == 1) {
                                target_cache[set_index].write_back[l] = 1;
                            }
                            break; 
                        }
                    }

                    // Hit 처리: 시간 갱신
                    if (hit) {
                        // ⭐ 가장 최근에 사용되었으므로 현재 시간으로 갱신
                        target_cache[set_index].lru_time[hit_way] = current_time;
                    }
                    // Miss 처리: LRU 교체
                    else {
                        if (is_icache) i_miss_count++;
                        else d_miss_count++;

                        int victim_way = -1;
                        unsigned long oldest_time = ULONG_MAX; // 최대값으로 초기화

                        // Victim 선정 로직
                        // 1. 빈 공간이 있으면 우선 사용
                        int found_empty = 0;
                        for (l = 0; l < assoc; l++) {
                            if (target_cache[set_index].valid[l] == 0) {
                                victim_way = l;
                                found_empty = 1;
                                break;
                            }
                        }

                        // 2. 빈 공간이 없으면 가장 오래된(LRU Time 최소) 블록 찾기
                        if (!found_empty) {
                            for (l = 0; l < assoc; l++) {
                                if (target_cache[set_index].lru_time[l] < oldest_time) {
                                    oldest_time = target_cache[set_index].lru_time[l];
                                    victim_way = l;
                                }
                            }
                        }

                        // [Write-Back Check]
                        if (target_cache[set_index].valid[victim_way] && target_cache[set_index].write_back[victim_way]) {
                            d_write_mem_count++;
                        }

                        // 새 블록 등록
                        target_cache[set_index].valid[victim_way] = 1;
                        target_cache[set_index].tag[victim_way] = tag_val;
                        
                        // 새 블록의 시간도 현재 시간으로 설정 (방금 들어왔으니까)
                        target_cache[set_index].lru_time[victim_way] = current_time;

                        // Dirty Bit 설정
                        if (!is_icache && cur_type == 1) {
                            target_cache[set_index].write_back[victim_way] = 1;
                        } else {
                            target_cache[set_index].write_back[victim_way] = 0;
                        }
                    }

                } /* Trace Loop End */

                /* ------------------------------------------------ */
                /* 3. 결과 저장 */
                /* ------------------------------------------------ */
                int row_i = assoc_idx;
                int row_d = assoc_idx + NUM_ASSOC;
                int col = (block_idx * NUM_CACHE) + cache_idx;

                miss[row_i][col] = (i_acc_count > 0) ? (double)i_miss_count / i_acc_count : 0.0;
                miss[row_d][col] = (d_acc_count > 0) ? (double)d_miss_count / d_acc_count : 0.0;

                writes[row_d][col] = d_write_mem_count;
                writes[row_i][col] = 0;

                i_totals[row_i][col] = i_acc_count;
                d_totals[row_d][col] = d_acc_count;

            } /* Cache Size Loop */
        } /* Block Size Loop */
    } /* Assoc Loop */
}

static void simulate_fifo(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {

    int i, k, l;
    int assoc_idx, block_idx, cache_idx;

    /* ------------------------------------------------------------------------- */
    /* 3중 루프: 모든 캐시 구성(Config)에 대해 시뮬레이션 수행 */
    /* ------------------------------------------------------------------------- */
    
    // Loop 1: Associativity (1, 2, 4, 8)
    for (assoc_idx = 0; assoc_idx < NUM_ASSOC; assoc_idx++) {
        int assoc = ASSOC_LIST[assoc_idx];

        // Loop 2: Block Size (8, 16, 32, 64, 128)
        for (block_idx = 0; block_idx < NUM_BLOCK; block_idx++) {
            int blk_size = BLOCK_SIZES[block_idx];

            // Loop 3: Cache Size (1KB ~ 16KB)
            for (cache_idx = 0; cache_idx < NUM_CACHE; cache_idx++) {
                int cache_size = CACHE_SIZES[cache_idx];

                /* ------------------------------------------------ */
                /* 1. 초기화 단계 (Initialization) */
                /* ------------------------------------------------ */
                
                // 캐시 메모리(Data/Tag) 0으로 초기화
                memset(icah_fifo, 0, sizeof(struct Block_FIFO) * MAX_SETS);
                memset(dcah_fifo, 0, sizeof(struct Block_FIFO) * MAX_SETS);

                // FIFO 포인터 0으로 초기화
                memset(icah_fifo_ptr, 0, sizeof(int) * MAX_SETS);
                memset(dcah_fifo_ptr, 0, sizeof(int) * MAX_SETS);

                // 통계 변수 초기화
                int i_miss_count = 0;
                int d_miss_count = 0;
                int i_acc_count = 0;
                int d_acc_count = 0;
                int d_write_mem_count = 0; // Dirty Block이 쫓겨날 때 메모리 쓰기

                /* ------------------------------------------------ */
                /* 2. 트레이스 시뮬레이션 (Trace Loop) */
                /* ------------------------------------------------ */
                for (k = 0; k < length; k++) {
                    
                    int cur_type = type[k];            // 0:Read, 1:Write, 2:Instruction
                    unsigned long cur_addr = addr[k];  // 주소

                    // 2-1. 주소 디코딩 (헬퍼 함수 사용)
                    unsigned long tag_val;
                    unsigned long set_index;
                    calculate_address_fields(cur_addr, blk_size, assoc, cache_size, 
                                             &tag_val, &set_index);

                    // 2-2. 타겟 캐시(I vs D) 설정
                    struct Block_FIFO *target_cache;
                    int *target_ptr;
                    int is_icache = (cur_type == 2);

                    if (is_icache) {
                        target_cache = icah_fifo;
                        target_ptr = icah_fifo_ptr;
                        i_acc_count++;
                    } else {
                        target_cache = dcah_fifo;
                        target_ptr = dcah_fifo_ptr;
                        d_acc_count++;
                    }

                    // 2-3. Hit 검사
                    int hit = 0;
                    int hit_way = -1;
                    
                    for (l = 0; l < assoc; l++) {
                        if (target_cache[set_index].valid[l] && target_cache[set_index].tag[l] == tag_val) {
                            hit = 1;
                            hit_way = l;
                            
                            // Hit & Write(D-Cache) -> Dirty Bit 설정
                            if (!is_icache && cur_type == 1) {
                                target_cache[set_index].write_back[l] = 1;
                            }
                            break; 
                        }
                    }

                    // 2-4. Miss 처리 (FIFO 정책 적용)
                    if (!hit) {
                        if (is_icache) i_miss_count++;
                        else d_miss_count++;

                        // Victim 선정: FIFO 포인터가 가리키는 곳
                        int victim_way = target_ptr[set_index];

                        // [Write-Back] 쫓겨나는 블록이 유효하고 Dirty라면 메모리 쓰기 발생
                        if (target_cache[set_index].valid[victim_way] && target_cache[set_index].write_back[victim_way]) {
                            d_write_mem_count++;
                        }

                        // 새 블록 등록 (Replace)
                        target_cache[set_index].valid[victim_way] = 1;
                        target_cache[set_index].tag[victim_way] = tag_val;

                        // 새 블록의 Dirty 상태 설정
                        // Write Miss면 가져와서 썼으므로 Dirty=1, 아니면 0
                        if (!is_icache && cur_type == 1) {
                            target_cache[set_index].write_back[victim_way] = 1;
                        } else {
                            target_cache[set_index].write_back[victim_way] = 0;
                        }

                        // FIFO 포인터 업데이트 (Circular Buffer)
                        target_ptr[set_index] = (victim_way + 1) % assoc;
                    }

                } /* Trace Loop End */

                /* ------------------------------------------------ */
                /* 3. 결과 저장 (Store Results) */
                /* ------------------------------------------------ */
                
                // 행(Row) 인덱스 계산
                int row_i = assoc_idx;
                int row_d = assoc_idx + NUM_ASSOC;

                // 열(Col) 인덱스 계산
                int col = (block_idx * NUM_CACHE) + cache_idx;

                // Miss Rate 저장 (0으로 나누기 방지 포함)
                miss[row_i][col] = (i_acc_count > 0) ? (double)i_miss_count / i_acc_count : 0.0;
                miss[row_d][col] = (d_acc_count > 0) ? (double)d_miss_count / d_acc_count : 0.0;

                // Write Count 저장 (D-Cache만 해당)
                writes[row_d][col] = d_write_mem_count;
                writes[row_i][col] = 0; 

                // BEST 정책 계산용 총 접근 수 저장
                i_totals[row_i][col] = i_acc_count;
                d_totals[row_d][col] = d_acc_count;

            } /* Cache Size Loop End */
        } /* Block Size Loop End */
    } /* Assoc Loop End */
}

static void simulate_new(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {


    /*MRU 사용*/
    int i, k, l;
    int assoc_idx, block_idx, cache_idx;
    
    /* 시간 추적을 위한 변수 (시뮬레이션 전체 시간) */
    unsigned long current_time = 0;

    /* 3중 루프: 모든 구성 시뮬레이션 */
    for (assoc_idx = 0; assoc_idx < NUM_ASSOC; assoc_idx++) {
        int assoc = ASSOC_LIST[assoc_idx];

        for (block_idx = 0; block_idx < NUM_BLOCK; block_idx++) {
            int blk_size = BLOCK_SIZES[block_idx];

            for (cache_idx = 0; cache_idx < NUM_CACHE; cache_idx++) {
                int cache_size = CACHE_SIZES[cache_idx];

                /* ------------------------------------------------ */
                /* 1. 초기화 (Initialization) */
                /* ------------------------------------------------ */
                
                // 캐시 메모리 및 포인터 초기화
                memset(icah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                memset(dcah_lru, 0, sizeof(struct Block_LRU) * MAX_SETS);
                
                // 통계 변수 초기화
                int i_miss_count = 0;
                int d_miss_count = 0;
                int i_acc_count = 0;
                int d_acc_count = 0;
                int d_write_mem_count = 0;

                /* ------------------------------------------------ */
                /* 2. 트레이스 시뮬레이션 (Trace Loop) */
                /* ------------------------------------------------ */
                for (k = 0; k < length; k++) {
                    current_time++; // ⭐ 시간 흐름 (매 접근마다 증가)

                    int cur_type = type[k];
                    unsigned long cur_addr = addr[k];

                    // 주소 디코딩 (헬퍼 함수)
                    unsigned long tag_val;
                    unsigned long set_index;
                    calculate_address_fields(cur_addr, blk_size, assoc, cache_size, 
                                             &tag_val, &set_index);

                    // 타겟 캐시 설정
                    struct Block_LRU *target_cache;
                    int is_icache = (cur_type == 2);

                    if (is_icache) {
                        target_cache = icah_lru;
                        i_acc_count++;
                    } else {
                        target_cache = dcah_lru;
                        d_acc_count++;
                    }

                    // Hit 검사
                    int hit = 0;
                    int hit_way = -1;

                    for (l = 0; l < assoc; l++) {
                        if (target_cache[set_index].valid[l] && target_cache[set_index].tag[l] == tag_val) {
                            hit = 1;
                            hit_way = l;
                            
                            // Hit & Write -> Dirty Bit 설정
                            if (!is_icache && cur_type == 1) {
                                target_cache[set_index].write_back[l] = 1;
                            }
                            break; 
                        }
                    }

                    // Hit 처리: 시간 갱신
                    if (hit) {
                        // ⭐ 가장 최근에 사용되었으므로 현재 시간으로 갱신
                        target_cache[set_index].lru_time[hit_way] = current_time;
                    }
                    // Miss 처리: LRU 교체
                    else {
                        if (is_icache) i_miss_count++;
                        else d_miss_count++;

                        int victim_way = -1;
                        unsigned long most_recent_time = 0;

                        // Victim 선정 로직
                        // 1. 빈 공간이 있으면 우선 사용
                        int found_empty = 0;
                        for (l = 0; l < assoc; l++) {
                            if (target_cache[set_index].valid[l] == 0) {
                                victim_way = l;
                                found_empty = 1;
                                break;
                            }
                        }

                        // 2. 빈 공간이 없으면 가장 오래된(LRU Time 최소) 블록 찾기
                        if (!found_empty) {
                            for (l = 0; l < assoc; l++) {
                                if (target_cache[set_index].lru_time[l] > most_recent_time) {
                                    most_recent_time = target_cache[set_index].lru_time[l];
                                    victim_way = l;
                                }
                            }
                        }

                        // [Write-Back Check]
                        if (target_cache[set_index].valid[victim_way] && target_cache[set_index].write_back[victim_way]) {
                            d_write_mem_count++;
                        }

                        // 새 블록 등록
                        target_cache[set_index].valid[victim_way] = 1;
                        target_cache[set_index].tag[victim_way] = tag_val;
                        
                        // 새 블록의 시간도 현재 시간으로 설정 (방금 들어왔으니까)
                        target_cache[set_index].lru_time[victim_way] = current_time;

                        // Dirty Bit 설정
                        if (!is_icache && cur_type == 1) {
                            target_cache[set_index].write_back[victim_way] = 1;
                        } else {
                            target_cache[set_index].write_back[victim_way] = 0;
                        }
                    }
                    // [추가 예시] 처음 10번의 접근에 대해서만 상태를 출력하여 검증
                    if (k < 10) {
                        print_new_cache_state(is_icache, set_index, assoc);
                    }

                } /* Trace Loop End */

                /* ------------------------------------------------ */
                /* 3. 결과 저장 */
                /* ------------------------------------------------ */
                int row_i = assoc_idx;
                int row_d = assoc_idx + NUM_ASSOC;
                int col = (block_idx * NUM_CACHE) + cache_idx;

                miss[row_i][col] = (i_acc_count > 0) ? (double)i_miss_count / i_acc_count : 0.0;
                miss[row_d][col] = (d_acc_count > 0) ? (double)d_miss_count / d_acc_count : 0.0;

                writes[row_d][col] = d_write_mem_count;
                writes[row_i][col] = 0;

                i_totals[row_i][col] = i_acc_count;
                d_totals[row_d][col] = d_acc_count;

            } /* Cache Size Loop */
        } /* Block Size Loop */
    } /* Assoc Loop */

}

static void print_results(const char *label,
                        const double miss[NUM_ROWS][NUM_COLS],
                        const int writes[NUM_ROWS][NUM_COLS]) {
    int i, j, k;

    printf("\nMissRate\n");
    for (i = 0; i < NUM_ROWS; i++) {

        if (i == 0) {
            printf("           ");
            for (k = 0; k < NUM_BLOCK; k++)
                printf("%s/%-4d                                 ", label, BLOCK_SIZES[k]);
            printf("\nI cache   ");
            for (k = 0; k < NUM_BLOCK; k++)
                for (j = 0; j < NUM_CACHE; j++) printf("%-7d", CACHE_SIZES[j]);
            printf("\n");
        }

        if (i == NUM_ASSOC) {
            printf("\n           ");
            for (k = 0; k < NUM_BLOCK; k++)
                printf("%s/%-4d                                 ", label, BLOCK_SIZES[k]);
            printf("\nD cache   ");
            for (k = 0; k < NUM_BLOCK; k++)
                for (j = 0; j < NUM_CACHE; j++) printf("%-7d", CACHE_SIZES[j]);
            printf("\n");
        }

        if (i % NUM_ASSOC == 0)      printf("Direct | ");
        else if (i % NUM_ASSOC == 1) printf("2  Way | ");
        else if (i % NUM_ASSOC == 2) printf("4  Way | ");
        else                         printf("8  Way | ");

        for (j = 0; j < NUM_COLS; j++)
            printf("%.4lf ", miss[i][j]);
        printf("\n");
    }

    printf("\nWrite Count\n");
    for (i = 0; i < NUM_ROWS; i++) {

        if (i == 0) {
            printf("           ");
            for (k = 0; k < NUM_BLOCK; k++)
                printf("%s/%-4d                          ", label, BLOCK_SIZES[k]);
            printf("\nI cache   ");
            for (k = 0; k < NUM_BLOCK; k++)
                for (j = 0; j < NUM_CACHE; j++)
                    printf("%6d", CACHE_SIZES[j]);
            printf("\n");
        }

        if (i == NUM_ASSOC) {
            printf("\n           ");
            for (k = 0; k < NUM_BLOCK; k++)
                printf("%s/%-4d                          ", label, BLOCK_SIZES[k]);
            printf("\nD cache   ");
            for (k = 0; k < NUM_BLOCK; k++)
                for (j = 0; j < NUM_CACHE; j++)
                    printf("%6d", CACHE_SIZES[j]);
            printf("\n");
        }

        if (i % NUM_ASSOC == 0)      printf("Direct | ");
        else if (i % NUM_ASSOC == 1) printf("2  Way | ");
        else if (i % NUM_ASSOC == 2) printf("4  Way | ");
        else                         printf("8  Way | ");

        for (j = 0; j < NUM_COLS; j++)
            printf("%5d ", writes[i][j]);

        printf("\n");
    }
}

static void print_best_results(
    const double lru_miss[NUM_ROWS][NUM_COLS], const int lru_writes[NUM_ROWS][NUM_COLS],
    const int lru_i_tot[NUM_ROWS][NUM_COLS], const int lru_d_tot[NUM_ROWS][NUM_COLS],
    const double fifo_miss[NUM_ROWS][NUM_COLS], const int fifo_writes[NUM_ROWS][NUM_COLS],
    const int fifo_i_tot[NUM_ROWS][NUM_COLS], const int fifo_d_tot[NUM_ROWS][NUM_COLS],
    int i_hit, int i_miss, int d_hit, int d_miss) {
    
    int cl, bl, as; // Loop variables: Cache Size, Block Size, Associativity

    // 1. 캐시 크기별로 최적의 설정을 찾습니다.
    for (cl = 0; cl < NUM_CACHE; cl++) {

        double min_cycles = DBL_MAX; // 최소 사이클 (초기값: 최대 실수)

        // 최적의 설정을 저장할 변수들
        const char *best_i_policy = "N/A";
        const char *best_d_policy = "N/A"; // I/D가 같은 정책을 쓰므로 사실 하나만 있어도 됨
        int best_block = 0;
        int best_assoc = 0;
        double best_i_missrate = 0.0;
        double best_d_missrate = 0.0;
        int best_d_writes = 0;
        double best_total_cycles = 0.0;

        // 2. 모든 블록 크기와 연관도 조합을 순회
        for (bl = 0; bl < NUM_BLOCK; bl++) {
            for (as = 0; as < NUM_ASSOC; as++) {
                
                // 인덱스 계산 (FIFO, LRU 공통)
                int row_i = as;
                int row_d = as + NUM_ASSOC;
                int col = (bl * NUM_CACHE) + cl;
                
                // -----------------------------------------------------
                // (A) LRU 성능 계산
                // -----------------------------------------------------
                double lru_i_mr = lru_miss[row_i][col];
                double lru_d_mr = lru_miss[row_d][col];
                int lru_i_total = lru_i_tot[row_i][col];
                int lru_d_total = lru_d_tot[row_d][col];

                double lru_i_misses = lru_i_total * lru_i_mr;
                double lru_i_hits   = lru_i_total - lru_i_misses;
                double lru_d_misses = lru_d_total * lru_d_mr;
                double lru_d_hits   = lru_d_total - lru_d_misses;

                double lru_cycles = (lru_i_hits * i_hit) + (lru_i_misses * i_miss) +
                                    (lru_d_hits * d_hit) + (lru_d_misses * d_miss);

                // LRU가 현재까지의 최소 사이클보다 작으면 갱신
                if (lru_cycles < min_cycles) {
                    min_cycles = lru_cycles;
                    best_i_policy = "LRU";
                    best_d_policy = "LRU";
                    best_block = BLOCK_SIZES[bl];
                    best_assoc = ASSOC_LIST[as];
                    best_i_missrate = lru_i_mr;
                    best_d_missrate = lru_d_mr;
                    best_d_writes = lru_writes[row_d][col];
                    best_total_cycles = lru_cycles;
                }

                // -----------------------------------------------------
                // (B) FIFO 성능 계산
                // -----------------------------------------------------
                double fifo_i_mr = fifo_miss[row_i][col];
                double fifo_d_mr = fifo_miss[row_d][col];
                int fifo_i_total = fifo_i_tot[row_i][col];
                int fifo_d_total = fifo_d_tot[row_d][col];

                double fifo_i_misses = fifo_i_total * fifo_i_mr;
                double fifo_i_hits   = fifo_i_total - fifo_i_misses;
                double fifo_d_misses = fifo_d_total * fifo_d_mr;
                double fifo_d_hits   = fifo_d_total - fifo_d_misses;

                double fifo_cycles = (fifo_i_hits * i_hit) + (fifo_i_misses * i_miss) +
                                     (fifo_d_hits * d_hit) + (fifo_d_misses * d_miss);

                // FIFO가 현재까지의 최소 사이클보다 작으면 갱신
                if (fifo_cycles < min_cycles) {
                    min_cycles = fifo_cycles;
                    best_i_policy = "FIFO";
                    best_d_policy = "FIFO";
                    best_block = BLOCK_SIZES[bl];
                    best_assoc = ASSOC_LIST[as];
                    best_i_missrate = fifo_i_mr;
                    best_d_missrate = fifo_d_mr;
                    best_d_writes = fifo_writes[row_d][col];
                    best_total_cycles = fifo_cycles;
                }
            }
        }

        // 3. 결과 출력 (Skeleton Code 형식 유지)
        printf("--- Cache Size: %d bytes ---\n", CACHE_SIZES[cl]);

        if (min_cycles == DBL_MAX) {
            printf("  I-Cache: No instruction accesses.\n");
            printf("  D-Cache: No data accesses.\n");
        } else {
            printf("  Best I-Cache: Policy=%-4s | Block=%-4d | Assoc=%-2d | MissRate=%.4f | Total Cycles=%.0f\n",
                   best_i_policy, best_block, best_assoc, best_i_missrate, best_total_cycles);

            printf("  Best D-Cache: Policy=%-4s | Block=%-4d | Assoc=%-2d | MissRate=%.4f | Writes=%-5d | Total Cycles=%.0f\n",
                   best_d_policy, best_block, best_assoc, best_d_missrate, best_d_writes, best_total_cycles);
        }

        printf("\n");
    }
}

static void print_two_cache_state(const char *policy, int is_icache, int index, int assoc) {
    printf("\n[Cache State Dump] Policy=%s | %s | index=%d | assoc=%d\n",
           policy,
           is_icache ? "I-Cache" : "D-Cache",
           index, assoc);

    // ---------------- LRU ----------------
    if (strcmp(policy, "LRU") == 0 || strcmp(policy, "lru") == 0) {

        struct Block_LRU *cache = is_icache ? icah_lru : dcah_lru;

        for (int i = 0; i < assoc; i++) {
            printf("  Way %-2d | valid=%d  tag=%lu  write_back=%d\n",
                   i,
                   cache[index].valid[i],
                   cache[index].tag[i],
                   cache[index].write_back[i]);
        }
    }

    // ---------------- FIFO ----------------
    else if (strcmp(policy, "FIFO") == 0 || strcmp(policy, "fifo") == 0) {

        struct Block_FIFO *cache = is_icache ? icah_fifo : dcah_fifo;
        int *ptr = is_icache ? icah_fifo_ptr : dcah_fifo_ptr;

        printf("  FIFO pointer = %d\n", ptr[index]);

        for (int i = 0; i < assoc; i++) {
            printf("  Way %-2d | valid=%d  tag=%lu  write_back=%d\n",
                   i,
                   cache[index].valid[i],
                   cache[index].tag[i],
                   cache[index].write_back[i]);
        }
    }

    printf("\n");
}

static void print_new_cache_state(int is_icache, int index, int assoc) {
    // 1. 헤더 출력
    printf("\n[NEW Policy State Dump] %s | index=%d | assoc=%d\n",
           is_icache ? "I-Cache" : "D-Cache",
           index, assoc);

    // 2. 캐시 배열 선택
    // MRU 구현 시 LRU 구조체와 배열(icah_lru, dcah_lru)을 재사용한다고 가정합니다.
    // 만약 simulate_new에서 별도의 배열(예: icah_new)을 선언해서 쓰셨다면 그걸로 바꿔주세요.
    struct Block_LRU *cache = is_icache ? icah_lru : dcah_lru;

    // 3. 각 Way의 상태 출력
    for (int i = 0; i < assoc; i++) {
        // 핵심: lru_time을 출력하여 교체 우선순위를 검증할 수 있게 함
        printf("  Way %-2d | valid=%d  tag=%lx  write_back=%d  lru_time=%lu\n",
               i,
               cache[index].valid[i],
               cache[index].tag[i],
               cache[index].write_back[i],
               cache[index].lru_time[i]); 
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    if (argc < 3 || (argc > 3 && argc < 7) || argc > 7)
        usage(argv[0]);

    int policy = 0;  
    int i_hit_c = 0, i_miss_c = 0, d_hit_c = 0, d_miss_c = 0;

    char *trace_file = NULL;

    if (!strcasecmp(argv[1], "FIFO")) {
        if (argc != 3) usage(argv[0]);
        policy = 1;
        trace_file = argv[2];
    }
    else if (!strcasecmp(argv[1], "LRU")) {
        if (argc != 3) usage(argv[0]);
        policy = 0;
        trace_file = argv[2];
    }
    else if (!strcasecmp(argv[1], "NEW")) {
        if (argc != 3) usage(argv[0]);
        policy = 3;
        trace_file = argv[2];
    }
    else if (!strcasecmp(argv[1], "BEST")) {
        if (argc != 7) usage(argv[0]);
        policy = 2;
        trace_file = argv[2];
        i_hit_c  = atoi(argv[3]);
        i_miss_c = atoi(argv[4]);
        d_hit_c  = atoi(argv[5]);
        d_miss_c = atoi(argv[6]);
    }
    else {
        usage(argv[0]);
    }

    int *type = NULL;
    unsigned long *addr = NULL;
    int length = 0;

    printf("Reading trace file: %s\n", trace_file);
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
    }

    else if (policy == 1) { 
        printf("Simulating FIFO policy...\n");
        double miss[NUM_ROWS][NUM_COLS] = {{0}};
        int writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int d_tot[NUM_ROWS][NUM_COLS]   = {{0}};

        simulate_fifo(type, addr, length, miss, writes, i_tot, d_tot);

        print_results("FIFO", miss, writes);
    }

    else if (policy == 3) { 
        printf("Simulating NEW policy...\n");
        double miss[NUM_ROWS][NUM_COLS] = {{0}};
        int writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int d_tot[NUM_ROWS][NUM_COLS]   = {{0}};

        simulate_new(type, addr, length, miss, writes, i_tot, d_tot);

        print_results("NEW", miss, writes);
    }

    else { 
        printf("Simulating LRU policy for BEST...\n");
        double lru_miss[NUM_ROWS][NUM_COLS] = {{0}};
        int lru_writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int lru_i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int lru_d_tot[NUM_ROWS][NUM_COLS]   = {{0}};

        double fifo_miss[NUM_ROWS][NUM_COLS] = {{0}};
        int fifo_writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int fifo_i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int fifo_d_tot[NUM_ROWS][NUM_COLS]   = {{0}};

        simulate_lru(type, addr, length, lru_miss, lru_writes, lru_i_tot, lru_d_tot);
        simulate_fifo(type, addr, length, fifo_miss, fifo_writes, fifo_i_tot, fifo_d_tot);

        printf("\n--- BEST Configuration Analysis ---\n");
        printf("Cycle Parameters: I(Hit/Miss) = %d/%d, D(Hit/Miss) = %d/%d\n\n",
               i_hit_c, i_miss_c, d_hit_c, d_miss_c);

        print_best_results(lru_miss, lru_writes, lru_i_tot, lru_d_tot,
                           fifo_miss, fifo_writes, fifo_i_tot, fifo_d_tot,
                           i_hit_c, i_miss_c, d_hit_c, d_miss_c);
    }

    free(type);
    free(addr);
    return 0;
}