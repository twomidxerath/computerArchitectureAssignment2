CC      = gcc
CFLAGS  = -O2 -Wall -lm
PROGRAM = cacheSim.out

# 'make'만 치면 실행되는 기본 타겟
all: $(PROGRAM)

# 실행 파일 생성 규칙
$(PROGRAM): main.c
	$(CC) $(CFLAGS) -o $(PROGRAM) main.c

# --- 테스트 케이스 정의 ---

test-trace1-fifo:
	@echo "Running trace1 FIFO..."
	@./$(PROGRAM) FIFO sample_input/trace1.txt > out_trace1_fifo.txt
	@diff -qwB out_trace1_fifo.txt sample_output/trace1-fifo.txt \
	    && echo "  PASS trace1-fifo" \
	    || (echo "  FAIL trace1-fifo"; exit 1)

test-trace1-lru:
	@echo "Running trace1 LRU..."
	@./$(PROGRAM) LRU sample_input/trace1.txt > out_trace1_lru.txt
	@diff -qwB out_trace1_lru.txt sample_output/trace1-lru.txt \
	    && echo "  PASS trace1-lru" \
	    || (echo "  FAIL trace1-lru"; exit 1)

test-trace1-best:
	@echo "Running trace1 BEST..."
	@./$(PROGRAM) BEST sample_input/trace1.txt 1 100 1 50 > out_trace1_best.txt
	@diff -qwB out_trace1_best.txt sample_output/trace1-best.txt \
	    && echo "  PASS trace1-best" \
	    || (echo "  FAIL trace1-best"; exit 1)

test-trace2-short-fifo:
	@echo "Running trace2-short FIFO..."
	@./$(PROGRAM) FIFO sample_input/trace2-short.txt > out_trace2_short_fifo.txt
	@diff -qwB out_trace2_short_fifo.txt sample_output/trace2-short-fifo.txt \
	    && echo "  PASS trace2-short-fifo" \
	    || (echo "  FAIL trace2-short-fifo"; exit 1)

test-trace2-short-lru:
	@echo "Running trace2-short LRU..."
	@./$(PROGRAM) LRU sample_input/trace2-short.txt > out_trace2_short_lru.txt
	@diff -qwB out_trace2_short_lru.txt sample_output/trace2-short-lru.txt \
	    && echo "  PASS trace2-short-lru" \
	    || (echo "  FAIL trace2-short-lru"; exit 1)

test-trace2-short-best:
	@echo "Running trace2-short BEST..."
	@./$(PROGRAM) BEST sample_input/trace2-short.txt 1 100 1 50 > out_trace2_short_best.txt
	@diff -qwB out_trace2_short_best.txt sample_output/trace2-short-best.txt \
	    && echo "  PASS trace2-short-best" \
	    || (echo "  FAIL trace2-short-best"; exit 1)

# 전체 테스트 실행 타겟
test: all \
      test-trace1-fifo \
      test-trace1-lru \
      test-trace1-best \
      test-trace2-short-fifo \
      test-trace2-short-lru \
      test-trace2-short-best
	@echo "========================================"
	@echo "        ALL TESTS PASSED"
	@echo "========================================"

clean:
	rm -f $(PROGRAM) out_*.txt