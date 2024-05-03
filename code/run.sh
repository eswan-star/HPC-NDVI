#!/usr/bin/env bash


#make sure we have environment set-up
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
export LD_LIBRARY_PATH=${GDAL_DIR}/lib64/:${TURBOJPEG_LIB}:${PROJ_DIR}/lib64:${OPENSSL_LIB}:${JSON_LIB}:${PAPI_LIB}:${TIFF_LIB}:${LD_LIBRARY_PATH}

parallel/parallel parallel/data/
