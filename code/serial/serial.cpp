
#include "Geotiff.hpp"

int main(int argc, char **argv)
{   
    //Get the filename "either a default or specified from command line"
    const char*pszFilename;
    if(argc==1)
        pszFilename = "NDVI_50_20231204_20231217_clip_20240407113518_2047798237.tif";
    else if(argc==2)
        pszFilename = argv[1];
    else
        return EINVAL;

    //Create a Geotiff object from
    Geotiff gtiff(pszFilename);

    //Print the desired metadata
    gtiff.GetRasterBand(1);

    //Get the desired data
    float* data; //ptr to float data of image
    GetArray2D(layerindex, data);

    //iterate over the data doing processing
    
    
    return 0;
}
