CC=icc

all: FDTD22 FDTD24

FDTD22: FDTD22_NOOPT.cpu FDTD22_SIMD.cpu FDTD22_SIMD_COLLAPSE.cpu FDTD22_SIMD_COLLAPSE_TILE.cpu

FDTD24: FDTD24_NOOPT.cpu FDTD24_SIMD.cpu FDTD24_SIMD_COLLAPSE.cpu FDTD24_SIMD_COLLAPSE_TILE.cpu

CFLAGS = -qopenmp -std=c99 -O3 -Wall

clean:
	rm -f *.cpu

FDTD22_NOOPT.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=2 -o FDTD22_NOOPT.cpu
FDTD22_SIMD.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=2 -DSIMD -o FDTD22_SIMD.cpu
FDTD22_SIMD_COLLAPSE.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=2 -DSIMD -DCOLLAPSE -o FDTD22_SIMD_COLLAPSE.cpu
FDTD22_SIMD_TILE.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=2 -DSIMD -DTILE -o FDTD22_SIMD_TILE.cpu
FDTD22_SIMD_COLLAPSE_TILE.cpu : main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=2 -DSIMD -DCOLLAPSE -DTILE -o FDTD22_SIMD_COLLAPSE_TILE.cpu


FDTD24_NOOPT.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=4 -o FDTD24_NOOPT.cpu
FDTD24_SIMD.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=4 -DSIMD -o FDTD24_SIMD.cpu
FDTD24_SIMD_COLLAPSE.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=4 -DSIMD -DCOLLAPSE -o FDTD24_SIMD_COLLAPSE.cpu
FDTD24_SIMD_TILE.cpu: main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=4 -DSIMD -DTILE -o FDTD24_SIMD_TILE.cpu
FDTD24_SIMD_COLLAPSE_TILE.cpu : main.c
	$(CC) main.c $(CFLAGS) -DSCHEME=4 -DSIMD -DCOLLAPSE -DTILE -o FDTD24_SIMD_COLLAPSE_TILE.cpu

