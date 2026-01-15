SHELL := /bin/bash
CFLAGS=-Wall -Wextra -O3 -ffast-math
LIBS=-lgsl -lgslcblas -lm
# NOTE: adjust `GSL_LIB_PATH` and `GSL_INC_PATH` according to your installation
# It is the path containing `libgsl.so.X` and `libgslcblas.so`
GSL_LIB_PATH=/usr/local/lib
GSL_INC_PATH=

BUILD_FOLDER=build

all: $(BUILD_FOLDER)/libbolt.so

$(BUILD_FOLDER)/libbolt.so: $(BUILD_FOLDER)/bolt.o
	gcc -shared -o $(BUILD_FOLDER)/libbolt.so $(BUILD_FOLDER)/bolt.o -I$(GSL_INC_PATH) -L$(GSL_LIB_PATH) -Wl,-rpath,$(GSL_LIB_PATH)/ $(LIBS)

$(BUILD_FOLDER)/bolt.o: bolt.c
	@if ! test -d $(BUILD_FOLDER); then \
		mkdir $(BUILD_FOLDER); \
	fi
	gcc $(CFLAGS) -c -fpic bolt.c -o $(BUILD_FOLDER)/bolt.o

.PHONY: test
test: bolt.c tests/matter_tk.c tests/matter_tk.py tests/distances.c tests/distances.py
	python3 ./tests/run_tests.py
	
bench: bolt.c tests/matter_tk
	@echo "Benchmarking..."
	gcc -o tests/matter_tk tests/matter_tk.c -lm -lgsl -lgslcblas -ggdb
	@time ./tests/matter_tk

clean:
	@if test -d $(BUILD_FOLDER); then \
		rm $(BUILD_FOLDER)/*; \
		rmdir $(BUILD_FOLDER); \
	fi