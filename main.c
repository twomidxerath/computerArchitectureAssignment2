#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define NUM_CACHE 5
#define NUM_BLOCK 5
#define NUM_ASSOC 4
#define NUM_ROWS 8
#define NUM_COLS 25

static const int CACHE_SIZES[NUM_CACHE] = {1024, 2048}; /* Modify this line */
static const int BLOCK_SIZES[NUM_BLOCK] = {64, 128};    /* Modify this line */
static const int ASSOC_LIST[NUM_ASSOC] = {1, 2};        /* Modify this line */


#define DO_NOT_USE_THIS 0               /* You should not use this */

/* LRU */
struct Block_LRU {
    char tag[DO_NOT_USE_THIS];          /* Modify this line */
    char valid[DO_NOT_USE_THIS];        /* Modify this line */
    char write_back[DO_NOT_USE_THIS];   /* Modify this line */
};
static struct Block_LRU icah_lru[1];    /* Modify this line */
static struct Block_LRU dcah_lru[1];    /* Modify this line */

/* FIFO */
struct Block_FIFO {
    char tag[DO_NOT_USE_THIS];          /* Modify this line */
    char valid[DO_NOT_USE_THIS];        /* Modify this line */
    char write_back[DO_NOT_USE_THIS];   /* Modify this line */
};
static struct Block_FIFO icah_fifo[1];  /* Modify this line */
static struct Block_FIFO dcah_fifo[1];  /* Modify this line */
static int icah_fifo_ptr[1];            /* Modify this line */
static int dcah_fifo_ptr[1];            /* Modify this line */


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

    /* ------------------------------------------------------------------------- */
    printf("Write your code here.\n");
    /* ------------------------------------------------------------------------- */

}

static void simulate_lru(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {

    /* ------------------------------------------------------------------------- */
    printf("Write your code here.\n");
    /* ------------------------------------------------------------------------- */

}

static void simulate_fifo(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {

    /* ------------------------------------------------------------------------- */
    printf("Write your code here.\n");
    /* ------------------------------------------------------------------------- */

}

static void simulate_new(int *type, unsigned long *addr, int length,
                        double miss[NUM_ROWS][NUM_COLS],
                        int writes[NUM_ROWS][NUM_COLS],
                        int i_totals[NUM_ROWS][NUM_COLS],
                        int d_totals[NUM_ROWS][NUM_COLS]) {

    /* ------------------------------------------------------------------------- */
    printf("Write your code here.\n");
    /* ------------------------------------------------------------------------- */

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
    
    int cl;

    for (cl = 0; cl < NUM_CACHE; cl++) {

        double best_i_time = DBL_MAX;
        double best_d_time = DBL_MAX;

        const char *best_i_policy = "N/A";
        const char *best_d_policy = "N/A";

        int best_i_block = 0, best_i_assoc = 0;
        int best_d_block = 0, best_d_assoc = 0;

        double best_i_missrate = 0.0;
        double best_d_missrate = 0.0;
        int best_d_writes = 0;

        /* ------------------------------------------------------------------------- */
        printf("Write your code here.\n");
        /* ------------------------------------------------------------------------- */

        printf("--- Cache Size: %d bytes ---\n", CACHE_SIZES[cl]);

        if (best_i_time == DBL_MAX)
            printf("  I-Cache: No instruction accesses.\n");
        else
            printf("  Best I-Cache: Policy=%-4s | Block=%-4d | Assoc=%-2d | MissRate=%.4f | Total Cycles=%.0f\n",
                   best_i_policy, best_i_block, best_i_assoc, best_i_missrate, best_i_time);

        if (best_d_time == DBL_MAX)
            printf("  D-Cache: No data accesses.\n");
        else
            printf("  Best D-Cache: Policy=%-4s | Block=%-4d | Assoc=%-2d | MissRate=%.4f | Writes=%-5d | Total Cycles=%.0f\n",
                   best_d_policy, best_d_block, best_d_assoc, best_d_missrate, best_d_writes, best_d_time);

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

    /* ------------------------------------------------------------------------- */
    printf("Write your code here.\n");
    /* ------------------------------------------------------------------------- */

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

        /* ------------------------------------------------------------------------- */
        printf("Write your code here.\n");
        /* ------------------------------------------------------------------------- */
        // simulate_lru( ... );
        
        print_results("LRU", miss, writes);
    }

    else if (policy == 1) { 
        printf("Simulating FIFO policy...\n");
        double miss[NUM_ROWS][NUM_COLS] = {{0}};
        int writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int d_tot[NUM_ROWS][NUM_COLS]   = {{0}};

        /* ------------------------------------------------------------------------- */
        printf("Write your code here.\n");
        /* ------------------------------------------------------------------------- */
        // simulate_fifo( ... );

        print_results("FIFO", miss, writes);
    }

    else if (policy == 3) { 
        printf("Simulating NEW policy...\n");
        double miss[NUM_ROWS][NUM_COLS] = {{0}};
        int writes[NUM_ROWS][NUM_COLS]  = {{0}};
        int i_tot[NUM_ROWS][NUM_COLS]   = {{0}};
        int d_tot[NUM_ROWS][NUM_COLS]   = {{0}};

        /* ------------------------------------------------------------------------- */
        printf("Write your code here.\n");
        /* ------------------------------------------------------------------------- */
        // simulate_new( ... );

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

        /* ------------------------------------------------------------------------- */
        printf("Write your code here.\n");
        /* ------------------------------------------------------------------------- */

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
