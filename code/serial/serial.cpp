#include "Geotiff.hpp"
#include <algorithm>
#include "papi.h"
#include <numeric>
#include <filesystem>
#include <regex>
#include <fstream>

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


void output_profile( long long int counters[4], long long int dt, double &MPT, double &IPB)
{
    long long total_cycles = counters[0];       // cpu cycles
    const long long total_instructions = counters[1]; // any
    const long long total_load_stores = counters[2];  // number of such instructions
    const long long total_l1d_misses = counters[3];   // number of access request to cache linereturn 0;

    const size_t flops = 1;
    //const size_t mem_ops = 1;
    double twall = static_cast<double>(dt) * 1.0e-9; // seconds
    const double IPC = static_cast<double>(total_instructions) / total_cycles;
    //const double OI =
    //               static_cast<double>(flops) / (total_load_stores * sizeof(float));
    //const double OI_theory =
    //                    static_cast<double>(flops) / (mem_ops * sizeof(float));
    const double float_perf = flops / twall * 1.0e-9; // Gflop/s
    MPT += static_cast<double>(total_load_stores)*sizeof(float) / twall; //mem op /s 
    IPB += static_cast<double>(total_instructions)/static_cast<double>(total_load_stores)/sizeof(float);

    std::cout<< "twall:                         " << twall << '\n';
    std::cout << "Total cycles:                 " << total_cycles << '\n';
    std::cout << "Total instructions:           " << total_instructions << '\n';
    std::cout << "Instructions per cycle (IPC): " << IPC << '\n';
    //std::cout << "Bandwidth: (Byte/s):          " << MPT << '\n';
    //std::cout << "Instructions per byte:        " << IPB << '\n';
    //std::cout << "Float performance:            " << float_perf <<'\n'; 
    return;
}

enum months{
    Jan,
    Feb,
    Mar,
    Apr,
    May,
    Jun,
    Jul,
    Aug,
    Sep,
    Oct,
    Nov,
    Dec
};


void write_to_file(map<string, float[12]> maxes, string filename)
{
    ofstream ofile(filename.c_str()); 
    ofile << "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec"<<endl;
    for( auto y : maxes){
        for(auto m : y.second)
            ofile << m <<" ";
        ofile<<endl;
    }
    ofile.close();
}


int main(int argc, char **argv)
{   
    //Get the filename "either a default or specified from command line"
    const char* pszFilename;
    const char* pszDir;
    if(argc==1)
        pszFilename = "NDVI_50_20231204_20231217_clip_20240407113518_2047798237.tif";
    else if(argc==2)
        pszDir = argv[1];
    else
        return EINVAL;

    //profiling stuff
    int event_set[2] = {PAPI_NULL, PAPI_NULL};
    int events[4] = {PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_LST_INS, PAPI_L1_DCM};
    long long int file_read_counters[4];
    long long int processing_counters[4];
    PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&event_set[0]);
    PAPI_add_events(event_set[0], events,4);
    PAPI_create_eventset(&event_set[1]);
    PAPI_add_events(event_set[1], events, 4); 
    long long int t0, t1;
    long long int times[4] = {0,0,0,0};

    //profile reading the files
    PAPI_start(event_set[0]);
    t0 = PAPI_get_real_nsec();
    map<string, float[12]> maxes;
    map<string, float[12]> mins;
    map<string, float[12]> means;
    map<string, float[12]> counts;

    //vector of vectors to store file data (not contiguous but avoids copy)
    vector<vector<float>> data;
    t0 = PAPI_get_real_nsec();
    for(const auto& pszFname : filesystem::recursive_directory_iterator(pszDir))
    {
        //get the file name
        string tmp = pszFname.path().string();
        pszFilename = tmp.c_str();
    
        //get a sense of what date this is from the file name
        regex pattern("202\\d{5}_");
        auto wbegin = sregex_iterator(tmp.begin(), tmp.end(), pattern);           
        auto wend = sregex_iterator();
        string year;
        int month;
        for(sregex_iterator i=wbegin; i!=wend; ++i)
        {
            smatch match = *i;
            string match_str = match.str();
            month = atoi(match_str.substr(4,2).c_str());
            string day = match_str.substr(6,2);
            year = match_str.substr(0,4);
            cout<<"month: "<<month<<", day:"<<day<<", year:"<<year<<endl;
        } 

        //If this year doesn't exist in map then add it
        if(maxes.find(year)==maxes.end())
        {
            memset(maxes[year],0,12);
            memset(means[year],0,12);
            memset(counts[year],0,12);
            memset(mins[year],0,12);
        }
       
        //Create object
        PAPI_start(event_set[0]);
        Geotiff gtiff(pszFilename);
        PAPI_stop(event_set[0], file_read_counters);
        data.push_back(gtiff.GetRasterBand(1));
    }   
    times[0] += PAPI_get_real_nsec() - t0;
    

    /////////do processing///////
    t0 = PAPI_get_real_nsec();
    //first need to convert counts to actual data values and throw out bounadry values
    size_t tot_size =  0;
    for(vector<float> vec: data)
    {
        
    }

    Image img(data, gtiff.GetDimensions());
    //iterate over the data doing processing
    float max_ndvi = -numeric_limits<float>::max();
    float min_ndvi = numeric_limits<float>::max();
    float mean = 0, count =0;
    PAPI_start(event_set[1]);
    t0 = PAPI_get_real_nsec();
    for(float &ndvi : img.vec){
        max_ndvi = ndvi>max_ndvi ? ndvi : max_ndvi;
        min_ndvi = ndvi<min_ndvi ? ndvi : min_ndvi;
        mean += ndvi;
        count+=1;
    }
    times[1] += PAPI_get_real_nsec() - t0;
    PAPI_stop(event_set[1], processing_counters);
    maxes[year][month-1] = max_ndvi>maxes[year][month-1] ? max_ndvi : maxes[year][month-1];        
    mins[year][month-1] = min_ndvi>mins[year][month-1] ? min_ndvi : mins[year][month-1];
    means[year][month-1] += mean;
    counts[year][month-1] += count;
    //std::cout<<"\n Max NDVI: "<<max_ndvi<<std::endl;
    }
    
    cout<<"\nNDVI trends:"<<endl;
    for(auto N : maxes){
        string year = N.first;
        cout<<"Year:"<<year<<"\n";
        for(int m=0;m<12;m++){
            if(counts[year][m]>0)
                cout<<"Month:"<<m+1<<", max ndvi: "<<N.second[m]<<", mean:"<<means[year][m]/counts[year][m]<<", min:"<<mins[year][m]<<endl;    
        }
    }
    
    double MPT=0, IPB=0; 
    std::cout<<"\nGDAL Geotiff initialization:"<<std::endl;
    output_profile(file_read_counters, times[0],MPT,IPB);
    
    //profile processing the data
    std::cout<<"\nData reduction:"<<std::endl;
    output_profile(processing_counters, times[1],MPT, IPB);
    
    double total_time = (static_cast<double>(times[0]) + static_cast<double>(times[1])) * 1e-9;
    cout<<"\nTotals:\n";
    cout << "Execution time:"<<total_time<<"s, GDAL init("<<times[0]/total_time*100*1e-9<<"%), Reduction("<<times[1]/total_time*100*1e-9<<"%)\n";     
    cout << "Traffic: (Byte/s):            " << MPT << '\n';
    cout << "Instructions per byte:        " << IPB << '\n';
 
    //std::cout<<"\nMax NDVI: "<<std::endl; 
    //for(float m : maxes) std::cout<<m<<", ";
    //std::cout<<std::endl;
    
    //Write out a text file which has timing data so python can use it
    write_to_file(maxes, "maxes.txt");
    write_to_file(means, "means.txt");
    write_to_file(mins, "mins.txt");

    return 1;
}
