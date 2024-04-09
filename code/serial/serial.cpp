#include "Geotiff.hpp"
#include <algorithm>
#include "papi.h"
#include <numeric>
#include <filesystem>
#include <regex>

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


void output_profile( long long int counters[4], long long int dt)
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

    std::cout<< "twall:                         " << twall << '\n';
    std::cout << "Total cycles:                 " << total_cycles << '\n';
    std::cout << "Total instructions:           " << total_instructions << '\n';
    std::cout << "Instructions per cycle (IPC): " << IPC << '\n';
    std::cout << "Float performance:            " << float_perf <<'\n';
   
    return;
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
    long long int times[2] = {0,0};

    //profile reading the files
    PAPI_start(event_set[0]);
    t0 = PAPI_get_real_nsec();
    map<string, float> monthly_ndvis;
    for(const auto& pszFname : filesystem::recursive_directory_iterator(pszDir))
    {
        std::cout<<pszFname<<std::endl;
        string tmp = pszFname.path().string();
        pszFilename = tmp.c_str();
    
        //get a sense of what date this is from the file name
        regex pattern("2024\\d{4}_");
        auto wbegin = sregex_iterator(tmp.begin(), tmp.end(), pattern);           
        auto wend = sregex_iterator();
        string month;
        for(sregex_iterator i=wbegin; i!=wend; ++i)
        {
            smatch match = *i;
            string match_str = match.str();
            month = match_str.substr(4,2);
            string day = match_str.substr(6,2);
            cout<<"month: "<<month<<", day:"<<day<<endl;
        }
        
        //Create object
        PAPI_start(event_set[0]);
        t0 = PAPI_get_real_nsec();
        Geotiff gtiff(pszFilename);
        vector<float> data = gtiff.GetRasterBand(1);
        times[0] += PAPI_get_real_nsec() - t0;
        PAPI_stop(event_set[0], file_read_counters);
        
        //do processing
        Image img(data, gtiff.GetDimensions());
        //iterate over the data doing processing
        float max_ndvi = numeric_limits<float>::min();
        PAPI_start(event_set[1]);
        t0 = PAPI_get_real_nsec();
        for(size_t i=0;i<img.vec.size();i++)
            max_ndvi = img.vec[i]>max_ndvi ? img.vec[i] : max_ndvi;
        monthly_ndvis[month]  = max_ndvi>monthly_ndvis[month] ? max_ndvi : monthly_ndvis[month];
        times[1] += PAPI_get_real_nsec() - t0;
        PAPI_stop(event_set[1], processing_counters);

        //std::cout<<"\n Max NDVI: "<<max_ndvi<<std::endl;
    }
    
    cout<<"\nMonthly max ndvi's:"<<endl;
    for(auto ndvi : monthly_ndvis)
        cout<<"Month: "<<ndvi.first<<", max ndvi:"<<ndvi.second<<endl;

    std::cout<<"\nReading files time:"<<std::endl;
    output_profile(file_read_counters, times[0]);
    
    //profile processing the data
    std::cout<<"\nData reduction time:"<<std::endl;
    output_profile(processing_counters, times[1]);
 
    //std::cout<<"\nMax NDVI: "<<std::endl; 
    //for(float m : maxes) std::cout<<m<<", ";
    //std::cout<<std::endl;
    
    return 1;
}
