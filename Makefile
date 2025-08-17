CFLAGS=-Wall -Wextra

libbolt.so: bolt.o
	gcc $(CFLAGS) -shared -o libbolt.so bolt.o -L/usr/local/lib -Wl,-rpath,/usr/local/lib/ -lgsl -lgslcblas -lm

bolt.o: bolt.c
	gcc -c -fpic bolt.c

clean:
	rm bolt.o
	rm libbolt.so