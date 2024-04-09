PROG=Geotiff
CC=x86_64-conda_cos7-linux-gnu-g++
CPPFLAGS= -g -Wall -I${CONDA_PREFIX}/include
LDFLAGS=-L${CONDA_PREFIX}/lib -lgdal -lm

all:
	@$(CC) -o serial serial.cpp  $(CPPFLAGS) $(LDFLAGS)

clean:
	@rm $(PROG)
