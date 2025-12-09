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

# TODO: make a test suite. For now, this compiles a minimal executable.
# In the future, I want to test the different features of Bolt.
.PHONY: test
test: bolt.c tests/hello_world.c tests/hello_world.py tests/distances.c tests/distances.py
	
	@echo "Compiling distances test..."
	gcc -o tests/distances tests/distances.c -lm -lgsl -lgslcblas -ggdb 
	@echo "Executing distances test..."
	@echo "Expected output:  z = 0.50, D_L = 2940.46 Mpc || z = 1.00, D_L = 6826.99 Mpc"
	@echo "Actual output C : `./tests/distances; echo Status: $$?`"
	@echo "Actual output Py: `python3 ./tests/distances.py; echo Status: $$?`"
	
	@echo "Compiling transfer function test..."
	gcc -o tests/hello_world tests/hello_world.c -lm -lgsl -lgslcblas -ggdb

	
	@echo "Executing transfer function test..."
	@echo "Expected output:  Hello, from Bolt! For h = 0.67, Omega_m = 0.319, \delta_m(k=0.1, z=0) = -65266.918234"
	@echo "Actual output C : `./tests/hello_world; echo Status: $$?`"
	@echo "Actual output Py: `python3 ./tests/hello_world.py; echo Status: $$?`"

bench: bolt.c tests/hello_world
	@echo "Benchmarking..."
	gcc -o tests/hello_world tests/hello_world.c -lm -lgsl -lgslcblas -ggdb
	@time ./tests/hello_world

clean:
	@if test -d $(BUILD_FOLDER); then \
		rm $(BUILD_FOLDER)/*; \
		rmdir $(BUILD_FOLDER); \
	fi