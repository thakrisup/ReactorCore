SRC = src/simulation_configuration.c src/main.c src/simulation_support.c
LFLAGS=-lm -lmpi
CFLAGS=-O3

.PHONY: archer2 local build

archer2: CC=cc
archer2: build

cirrus: CC=gcc
cirrus: build

build: 
	$(CC) -g -o reactor $(SRC) $(CFLAGS) $(LFLAGS)

