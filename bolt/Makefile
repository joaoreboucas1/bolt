CFLAGS=-Wall -Wextra
LIBS=-lgsl -lgslcblas -lm
# NOTE: adjust `GSL_LIB_PATH` according to your installation
# It is the path containing `libgsl.so.X` and `libgslcblas.so`
GSL_LIB_PATH=/usr/local/lib

BUILD_FOLDER=build

$(BUILD_FOLDER)/libbolt.so: $(BUILD_FOLDER)/bolt.o
	gcc $(CFLAGS) -shared -o $(BUILD_FOLDER)/libbolt.so $(BUILD_FOLDER)/bolt.o -L$(GSL_LIB_PATH) -Wl,-rpath,$(GSL_LIB_PATH)/ $(LIBS)

$(BUILD_FOLDER)/bolt.o: bolt.c
	@if ! test -d $(BUILD_FOLDER); then \
		mkdir $(BUILD_FOLDER); \
	fi
	gcc -c -fpic bolt.c -o $(BUILD_FOLDER)/bolt.o

clean:
	@if test -d $(BUILD_FOLDER); then \
		rm $(BUILD_FOLDER)/*; \
		rmdir $(BUILD_FOLDER); \
	fi