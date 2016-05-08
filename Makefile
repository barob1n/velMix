


all: velMix.exe

velMix.exe: velMix.o segyIO_class.o
	gcc -O2 -o  velMix.exe velMix.o segyIO_class.o -lm -fopenmp

velMix.o: velMix.c
	gcc -O2 -c velMix.c -lm -fopenmp

segyIO_class.o: segyIO_class.c
	gcc -O2 -c segyIO_class.c


