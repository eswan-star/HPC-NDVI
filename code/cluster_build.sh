#Load needed modules
spack load python@3.11.6%gcc@7.3.1 cmake@3.27.7%gcc@7.3.1 libtiff@4.5.1 sqlite@3.43.2%gcc@7.3.1 json-c@0.16 hdf5/43gozkm zlib-ng@2.1.4%gcc@7.3.1

export LD_LIBRARY_PATH=/shared/software/spack/opt/spack/linux-amzn2-skylake_avx512/gcc-7.3.1/libtiff-4.5.1-fiyvp24cb4yzjg5cw6bkvh2x4l2s4ehy/lib64/:${LD_LIBRARY_PATH}

TOP_DIR=`pwd`
PROJ_DIR=${TOP_DIR}/proj-9.3.1
GDAL_DIR=${TOP_DIR}/gdal-3.8.5

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
fi

# return to containg directory
cd ${TOP_DIR}

#now check for gdal and if build if neededd
if [ ! -d ${GDAL_DIR} ]
then
	wget https://github.com/OSGeo/gdal/releases/download/v3.8.5/gdal-3.8.5.tar.gz
	tar -xzvf gdal-3.8.5.tar.gz
	rm gdal-3.8.5.tar.gz
	mkdir ${GDAL_DIR}/build -p
fi

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
