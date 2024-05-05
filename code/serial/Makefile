CXX=g++
#PAPI_DIR ?=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/papi-6.0.0.1-jevzvgsbwnpxam5t7obqxl45ztwxwvhi
EIGEN_INC=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/eigen-3.3.7-pwe4bmtruhdeqccai3qasfyctdmb3qlj/include
JSON_LIB=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/json-c-0.16-2g6i65hxqveht3gwzu3qymsq3js6jp4j/lib64
OPENSSL_LIB=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/openssl-3.1.3-lbbwfmbwca7cmvxrlsdxwq7ax5zgtch4/lib64
CPPFLAGS=-g -Wall -I${GDAL_DIR}/include -I${PAPI_DIR}/include -I${EIGEN_INC} 
LDFLAGS=-L${PROJ_DIR}/lib64 -L${GDAL_DIR}/lib64 -L${PAPI_DIR}/lib -L${JSON_LIB} -L${OPENSSL_LIB} -L${TURBOJPEG_LIB}
LFLAGS=-lgdal -lm -lpapi -lstdc++fs -lturbojpeg

all: serial

serial:  
	$(CXX) -o $@ serial.cpp ${CPPFLAGS} ${LDFLAGS} ${LFLAGS}

clean:
	rm $(PROG)
