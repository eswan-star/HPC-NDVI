PROG=serial
CC=x86_64-conda_cos7-linux-gnu-g++
PAPI_DIR=${HOME}/team18_2024/code/papi-7.1.0
CPPFLAGS=-g -Wall -I${CONDA_PREFIX}/include -I${PAPI_DIR}/include
LDFLAGS=-L${CONDA_PREFIX}/lib -L${PAPI_DIR}/lib 
LFLAGS=-lgdal -lm -lpapi

all:
	@$(CC) serial.cpp  $(CPPFLAGS) $(LDFLAGS) -o $(PROG) ${LFLAGS}

clean:
	@rm $(PROG)
