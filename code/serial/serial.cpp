
#include "Geotiff.hpp"
#include <algorithm>

//Just a wrapper to index 2D image stored as flat array
class Image
{
    private:
        int Nrow;
        int Ncol;
       
    public:
        vector<float> vec;
        float *data;
        
        Image(vector<float> &_data, int _dims[3])
        {
            //get image dimensions
            Nrow = _dims[0];
            Ncol = _dims[1];
            data = _data.data();    
            vec = _data;
        }

        float operator()(int row, int col)
        {
            return data[row*Ncol+col]; //row major order
        }


};


void func()
{
    
    return;
}


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
    int* dims = gtiff.GetDimensions();
    for(int d=0;d<3;d++)
        std::cout<<dims[d]<<", ";
    std::cout<<std::endl;

    //Get the desired data
    int layerindex=0;
    vector<float> data = gtiff.GetRasterBand(layerindex);
    Image img(data, dims);

    //iterate over the data doing processing
    float max_ndvi = numeric_limits<float>::max();
    for(size_t i=0;i<img.vec.size();i++)
        max_ndvi = img.vec[i]>max_ndvi ? img.vec[i] : max_ndvi;
    
    std::cout<<"Max NDVI: "<<max_ndvi<<std::endl;
    
    return 0;
}
