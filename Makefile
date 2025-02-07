CXX = g++
CXXFLAGS = -Wall -std=c++17 -fopenmp -DLIBSAIS_OPENMP -Iinclude -Iext/libsais
CC = g++
CFLAGS = -fopenmp -DLIBSAIS_OPENMP -Iinclude -Iext/libsais

SRCS = main.cpp src/utils.cpp src/mem_finder.cpp src/sequence_split_align.cpp ext/SW/ssw.cpp ext/SW/ssw_cpp.cpp

ifdef DEBUG
	CXXFLAGS += -O0 -g -DDEBUG
	CFLAGS += -O0 -g
else
	CXXFLAGS += -O3
	CFLAGS += -O3
endif

ifeq ($(OS),Windows_NT)
else
	CXXFLAGS += -lpthread -pthread -lrt
	SRCS += src/thread_pool.cpp src/thread_condition.cpp
endif

ifdef M64
	CXXFLAGS += -DM64
	CFLAGS += -DM64
	LIBSAIS_OBJS = ext/libsais/libsais64.o
	LIBSAIS_SRC = ext/libsais/libsais64.c
	LIBSAIS_HDR = ext/libsais/libsais64.h
else
	LIBSAIS_OBJS = ext/libsais/libsais.o
	LIBSAIS_SRC = ext/libsais/libsais.c
	LIBSAIS_HDR = ext/libsais/libsais.h
endif

OBJS = $(SRCS:.cpp=.o) $(LIBSAIS_OBJS)

FMAlign2: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o FMAlign2

utils.o: src/utils.cpp include/utils.h include/common.h include/kseq.h
	$(CXX) $(CXXFLAGS) -c src/utils.cpp -o $@

mem_finder.o: src/mem_finder.cpp include/common.h $(LIBSAIS_HDR) include/thread_pool.h
	$(CXX) $(CXXFLAGS) -c src/mem_finder.cpp -o $@

sequence_split_align.o: src/sequence_split_align.cpp include/common.h
	$(CXX) $(CXXFLAGS) -c src/sequence_split_align.cpp -o $@

$(LIBSAIS_OBJS): $(LIBSAIS_SRC) $(LIBSAIS_HDR)
	$(CC) $(CFLAGS) -c $(LIBSAIS_SRC) -o $@

ifeq ($(OS),Windows_NT)
else
thread_pool.o: src/thread_pool.cpp include/thread_pool.h include/thread_condition.h
	$(CXX) $(CXXFLAGS) -c src/thread_pool.cpp -o $@

thread_condition.o: src/thread_condition.cpp include/thread_condition.h
	$(CXX) $(CXXFLAGS) -c src/thread_condition.cpp -o $@
endif

ssw.o: ext/SW/ssw.cpp ext/SW/ssw.h
	$(CXX) $(CXXFLAGS) -c ext/SW/ssw.cpp -o $@

ssw_cpp.o: ext/SW/ssw_cpp.cpp ext/SW/ssw_cpp.h ext/SW/ssw.h
	$(CXX) $(CXXFLAGS) -c ext/SW/ssw_cpp.cpp -o $@
	
main.o: main.cpp include/utils.h include/common.h include/mem_finder.h include/sequence_split_align.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o $@

clean:
ifeq ($(OS),Windows_NT)
	del /F /Q $(OBJS) FMAlign2.exe 2> NUL || cmd /c exit 0
else
	rm -f $(OBJS) FMAlign2
endif
