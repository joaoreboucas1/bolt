libbolt.so: bolt.o
	gcc -shared -o libbolt.so bolt.o

bolt.o: bolt.c
	gcc -c -fpic bolt.c

clean:
	rm bolt.o
	rm libbolt.so