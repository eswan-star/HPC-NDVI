#!/usr/bin/env bash

#make sure we have spack set-up
. /shared/software/spack/share/spack/setup-env.sh

#Load needed modules
spack load python@3.11.6%gcc@7.3.1 cmake@3.27.7%gcc@7.3.1 libtiff@4.5.1 sqlite@3.43.2%gcc@7.3.1 json-c@0.16 hdf5/43gozkm 
spack load zlib-ng@2.1.4%gcc@7.3.1 papi@6.0.0.1 eigen@3.3.7 openssl@3.1.3%gcc@7.3.1 libjpeg-turbo@3.0.0%gcc@7.3.1 

#get whereabouts
TOP_DIR=`pwd`
export PROJ_DIR=${TOP_DIR}/proj-9.3.1
export GDAL_DIR=${TOP_DIR}/gdal-3.8.5

#point environment variables to installed packages
export PAPI_DIR=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/papi-6.0.0.1-jevzvgsbwnpxam5t7obqxl45ztwxwvhi
export PAPI_LIB=${PAPI_DIR}/lib/
export PAPI_BIN=${PAPI_DIR}/bin/
export TIFF_LIB=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/libtiff-4.5.1-fiyvp24cb4yzjg5cw6bkvh2x4l2s4ehy/lib64
export JSON_LIB=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/json-c-0.16-2g6i65hxqveht3gwzu3qymsq3js6jp4j/lib64
export OPENSSL_LIB=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/openssl-3.1.3-lbbwfmbwca7cmvxrlsdxwq7ax5zgtch4/lib64
export TURBOJPEG_LIB=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/libjpeg-turbo-3.0.0-3bwcg2epbokfwg23jqkvcsmq6qntntkl/lib64
export LD_LIBRARY_PATH=${TURBOJPEG_LIB}:${PROJ_DIR}/lib64:${OPENSSL_LIB}:${JSON_LIB}:${PAPI_LIB}:${TIFF_LIB}:${LD_LIBRARY_PATH}

export CMAKE_BIN=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/cmake-3.27.7-i7l5dirfl6fbqdkctneygfpqwgh27rzy/bin
export EIGEN_INC=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/eigen-3.3.7-pwe4bmtruhdeqccai3qasfyctdmb3qlj/include
export PATH=:${PAPI_BIN}:${CMAKE_BIN}:${PATH}

#install proj if needed
if [ ! -d ${PROJ_DIR} ] 
then
	wget https://download.osgeo.org/proj/proj-9.3.1.tar.gz
	tar -xzvf proj-9.3.1.tar.gz
	rm proj-9.3.1.tar.gz
	mkdir ${PROJ_DIR}/build -p
	cd ${PROJ_DIR}/build
	cmake .. -DCMAKE_INSTALL_PREFIX=${PROJ_DIR} -DCMAKE_BUILD_TYPE=Release
	cmake --build .
	cmake --build . --target install
	
	# return to containg directory
	cd ${TOP_DIR}
fi

#now check for gdal and if build if needed
if [ ! -d ${GDAL_DIR} ]
then
	wget https://github.com/OSGeo/gdal/releases/download/v3.8.5/gdal-3.8.5.tar.gz
	tar -xzvf gdal-3.8.5.tar.gz
	rm gdal-3.8.5.tar.gz
	mkdir ${GDAL_DIR}/build -p

	#now build and install gdal
	cd ${GDAL_DIR}/build

	if [ -f CMakeCache.txt ]
	then
		rm CMakeCache.txt
	fi

	cmake .. -DCMAKE_INSTALL_PREFIX=${GDAL_DIR} -DGDAL_BUILD_OPTIONAL_DRIVERS=OFF -DOGR_BUILD_OPTIONAL_DRIVERS=OFF\
					-DPROJ_DIR=${PROJ_DIR} -DPROJ_INCLUDE_DIR=${PROJ_DIR}/include -DPROJ_LIBRARY=${PROJ_DIR}/lib64/libproj.so  -DCMAKE_BUILD_TYPE=Release\
					-DENABLE_DEFLATE64=OFF
	cmake --build . 
	cmake --build . --target install

	#go back where we started
	cd ${TOP_DIR}
fi

#build targets
cd parallel/
make

cd ${TOP_DIR}

