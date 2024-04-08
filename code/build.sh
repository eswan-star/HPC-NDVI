#wget https://download.osgeo.org/proj/proj-9.3.1.tar.gz
#wget https://github.com/OSGeo/gdal/releases/download/v3.8.5/gdal-3.8.5.tar.gz
#tar -xzvf gdal-3.8.5.tar.gz
#tar -xzvf proj-9.3.1.tar.gz

#Load needed modules
spack load python@3.11.6%gcc@7.3.1 cmake@3.27.7%gcc@7.3.1 libtiff@4.5.1 sqlite@3.43.2%gcc@7.3.1 json-c@0.16 hdf5/43gozkm zlib-ng@2.1.4%gcc@7.3.1

CUR_DIR=`pwd`
PROJ_DIR=${HOME}/team18_2024/code/proj-9.3.1
GDAL_DIR=${HOME}/team18_2024/code/gdal-3.8.5

#first need to install proj
#cd ${PROJ_DIR}/build
#cmake .. -DCMAKE_INSTALL_PREFIX=${PROJ_DIR} -DCMAKE_BUILD_TYPE=Release
#cmake --build .
#cmake --build . --target install

#now build and install gdal
cd ${GDAL_DIR}/build
rm CMakeCache.txt
cmake .. -DCMAKE_INSTALL_PREFIX=${GDAL_DIR} -DGDAL_BUILD_OPTIONAL_DRIVERS=OFF -DOGR_BUILD_OPTIONAL_DRIVERS=OFF\
                -DPROJ_DIR=${PROJ_DIR} -DPROJ_INCLUDE_DIR=${PROJ_DIR}/include -DPROJ_LIBRARY=${PROJ_DIR}/lib64 -DCMAKE_BUILD_TYPE=Release
cmake --build . 
cmake --build . --target install

#go back where we started
cd ${CUR_DIR}



