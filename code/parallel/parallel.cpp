#include <algorithm>
#include <papi.h>
#include <numeric>
#include <experimental/filesystem>
#include <regex>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <x86intrin.h>
#include <stdalign.h>
#include <cstring>
#include "Geotiff.hpp"
#include "transforms.hpp"

using namespace std;

//Just a wrapper to index 2D image stored as flat array
class Image
{
    public:
        vector<float> vec;
        float *data;
        float mean;
	int Nrow;
	int Ncol;
	
	Image(void)
	{
	    data = NULL;
	    mean = numeric_limits<float>::min();
	    Nrow=-1;
	    Ncol=-1;
	}

        Image(vector<float> &_data, int _dims[2])
        {
            //get image dimensions
            Nrow = _dims[0];
            Ncol = _dims[1];
            data = _data.data();    
            vec = _data;
	    mean = -numeric_limits<float>::max();
        }

	Image(vector<float> &_data, int Nr, int Nc)
	{
            Nrow = Nr;
            Ncol = Nc;
            data = _data.data();    
            vec = _data;
	    mean = -numeric_limits<float>::max();
	}
	
	void convert_to_delta(void)
	{
	    int cnt = 0;
            mean = 0;
	    for(int i=0; i<Nrow;++i){
		for(int j=0;j<Ncol;++j){
		    int idx = i*Ncol+j;
	            if(data[idx]<255){
                        data[idx] = (data[idx]-125)/125.0;
                        mean += data[idx]; 
			++cnt;	
                    }
	            else
			data[idx] = 0;
		}
	    }
	    mean /= cnt;
	    for(int i=0;i<Nrow;++i)
	        for(int j=0;j<Ncol;++j){
		    int idx = i*Ncol+j;
		    data[idx] = (data[idx]-mean)/mean;
		}
	    return;
	}

        float &operator()(int row, int col)
        {
            return data[row*Ncol+col]; 
        }

	int inline size(void)
	{
	    return Nrow*Ncol; 
	}
};


void write_binary(float* data, int Nt, int Nr, int Nc, string fname)
{
    ofstream file(fname.c_str(), ios::out|ios::binary);
    file.write((char*)&Nt, sizeof(int));
    file.write((char*)&Nr, sizeof(int));
    file.write((char*)&Nc, sizeof(int));
    file.write((char*)data, Nr*Nc*Nt*sizeof(float));
    file.close();
}


//just organizes images by date
struct  DateContainer
{
   map<int, map<int, map<int, Image>>> images;
   size_t size = 0;

   //accesses image(s) for particular year/month/day
   Image &operator()(int year, int month, int day){
        return images[year][month][day];
   }

   //accesses image(s) for particular year/month/day
   Image &operator()(string year, string month, string day){
        return images[atoi(year.c_str())][atoi(month.c_str())-1][atoi(day.c_str())];
   }

   //add image while keeping count of how many you've added
   void add(string year, string month, string day, Image & addv)
   {
        images[atoi(year.c_str())][atoi(month.c_str())-1][atoi(day.c_str())] = addv;
        size++;
        return;    
   }   
};

void write_analysis(long long int Pk_counters[5], float twall, float twall_rel[4],  float traffic, 
		float Tp_sec, float Ts_sec, size_t Ws, float Ts, size_t Wp, float Tp, String filename)
{
    ofstream ofile(filename.c_str()); 
    ofile << "Wall time(s):            "<< twall << endl;
    //Gflop/s
    float gflops = (Pk_counters[0] + Pk_counters[4]) / twall;
    ofile << "Performance (Gflops/s):  "<<gflops<<endl;
    //percent of peak
    ofile << "frac peak perf:          "<<gflops/3340<<"%"<<endl;
    //Traffic
    ofile<< "Mem traffic (GB/s):       "<< traffic <<endl;
    //p strong scaling
    ofile << "p:                       "<<  Tp_sec/Ts_sec << endl;
    ofile << "Sp_strong:               "<< Tp/Ts<<endl;
    ofile << "Ws,Ts:                   "<< Ws<<","<<Ts<<endl;
    ofile << "Wp,Tp:                   "<< Wp<<","<<Tp<<endl;

    ofile.close();
   
}

void output_profile( long long int counters[5], long long int dt)
{
    float sflops = static_cast<float>(counters[0])*1e-9;       
    float vflops = static_cast<float>(counters[5])*1e-9;
    long long int traffic = counters[1]*sizeof(float); 
    long long int total_l1d_cache_misses = counters[2]; 
    
    double twall = static_cast<double>(dt) * 1.0e-9; // seconds
    cout << "Wall time(s):                 " << twall << '\n';
    cout << "SFlops:                        " << sflops << '\n';
    cout << "VFlops:                       " << vflops << '\n';
    cout << "Traffic (bytes):              " << traffic << '\n';
    cout << "Operational Intensity:        " << static_cast<double>(sflops)/traffic << '\n';
    cout << "L1d Cache Misses:             " << total_l1d_cache_misses << '\n'; 
    return;
}


int main(int argc, char **argv)
{   
    //Get the filename "either a default or specified from command line"
    const char* pszFilename;
    const char* pszDir;

    if(argc==2)
        pszDir = argv[1];
    else
        return EINVAL;

    //profiling stuff
    auto ret  = PAPI_library_init(PAPI_VER_CURRENT);
    if( ret != PAPI_VER_CURRENT){
        if(ret>0)
		throw::runtime_error("PAPI library mismatch!\n");
	else{
	    cout<<"PAPI_EINVAL?:  "<<(ret==PAPI_EINVAL)<<endl;
	    cout<<"PAPI_ENOMEM?:  "<<(ret==PAPI_ENOMEM)<<endl;
	    cout<<"PAPI_ESBSTR?:  "<<(ret==PAPI_ESBSTR)<<endl;
	    cout<<"PAPI_ECMP?:    "<<(ret==PAPI_ECMP)<<endl;
	    cout<<"PAPI ESYS?:    "<<(ret==PAPI_ESYS)<<endl;
	    throw::runtime_error("Problem initializing PAPI\n");
	}
    }
    int event_set = PAPI_NULL;
    int events[5]    = {PAPI_SP_OPS, PAPI_LST_INS, PAPI_L1_DCM, PAPI_L1_ICM, PAPI_VEC_SP};
    long long int file_read_counters[5]  = {0,0,0,0,0};
    long long int processing_counters[5] = {0,0,0,0,0};
    long long int Pk_counters[5]         = {0,0,0,0,0};
    long long int Pk_counters_serial[5]  = {0,0,0,0,0};
    long long int output_counters[5]     = {0,0,0,0,0};
    PAPI_create_eventset(&event_set);
    PAPI_add_events(event_set, events, 5);
    long long int t0;
    long long int times[5] = {0,0,0,0,0};
    int lst_Nr, lst_Nc;

    /********************************read the files*************************/
    cout<<"Read the files"<<endl;
    DateContainer data;
    size_t tot_nel_read=0, Nt=0;
    t0 = PAPI_get_real_nsec();
    PAPI_start(event_set);
    for(const auto& pszFname : experimental::filesystem::recursive_directory_iterator(pszDir))
    {
        //get the file name
        string tmp = pszFname.path().string();
        pszFilename = tmp.c_str();
    
        //get a sense of what date this is from the file name
        regex pattern("20\\d{6}_");
        auto wbegin = sregex_iterator(tmp.begin(), tmp.end(), pattern);           
        auto wend = sregex_iterator();
        string year, month, day;
        for(sregex_iterator i=wbegin; i!=wend; ++i)
        {
            smatch match = *i;
            string match_str = match.str();
            year = match_str.substr(0,4); 
            month = match_str.substr(4,2);
            day = match_str.substr(6,2);
            //cout<<"month: "<<month<<", day:"<<day<<", year:"<<year<<endl;
        } 
       
        //Create object
        Geotiff gtiff(pszFilename);
	Image img(gtiff.GetRasterBand(1), gtiff.GetDimensions());
	img.convert_to_delta();
        lst_Nr = img.Nrow;
	lst_Nc = img.Ncol;
	data.add(year, month, day, img); 
        tot_nel_read += data(year, month, day).size();
        Nt++;
    }   
    PAPI_stop(event_set, file_read_counters);
    //PAPI_accum(event_set, file_read_counters);
    times[0] += PAPI_get_real_nsec() - t0;
    
    cout<<"starting processing: N:"<<data.size<<endl;
    /***************************do pre-processing******************************/
    //get the date/loop over images
    t0 = PAPI_get_real_nsec();
    PAPI_start(event_set
    alignas(64) float Pk[lst_Nr*lst_Nc*Nt]; //(float*)malloc(lst_Nr*lst_Nc*Nt*sizeof(float));
    vector<tuple<int,int,int>> dates;
    size_t nt=0;
    for(auto &yearv : data.images){
        int year = yearv.first;
        for(auto &monthv : yearv.second){
            int month = monthv.first;
            for(auto &dayv : monthv.second){
                int day = dayv.first;
                //for(float &val : dayv.second){  
                //}
		if(lst_Nr!=dayv.second.Nrow or lst_Nc!=dayv.second.Ncol){
		    free(Pk);
		    throw::runtime_error("Images inconsistent-sizes across time!\n");
		}
		lst_Nr = dayv.second.Nrow;
		lst_Nc = dayv.second.Ncol;
                memcpy(&Pk[lst_Nr*lst_Nc*nt], dayv.second.data, lst_Nr*lst_Nc*sizeof(float));
		dates.push_back(make_tuple(year,month,day));
                nt++;
            }
        }
    }
    times[1] += PAPI_get_real_nsec() - t0;
    //PAPI_accum(event_set, processing_counters);
    PAPI_stop(event_set, processing_counters);

    /*********************************now do the correlation function analysis*****************************/
    cout<<"Do analysis Nt:"<<Nt<<", Nr:"<<lst_Nr<<", Nc:"<<lst_Nc<<endl;
    t0 = PAPI_get_real_nsec();
    PAPI_start(event_set);
    get_Pk_kernel(Pk, Nt, lst_Nr, lst_Nc);
    times[2] += PAPI_get_real_nsec() - t0;
    PAPI_stop(event_set, Pk_counters); 

    //now do serial
    t0 = PAPI_get_real_nsec();
    PAPI_start(event_set);
    get_Pk_kernel_serial(Pk, Nt, lst_Nr, lst_Nc);
    times[3] += PAPI_get_real_nsec() - t0;
    PAPI_stop(event_set, Pk_counters_serial); 

    /*********************************now do the output**************************/
    t0 = PAPI_get_real_nsec();
    PAPI_start(event_set);
    //Write out a text file which has timing data so python can use it 
    write_binary(Pk, Nt, lst_Nr, lst_Nc, "Pk.bin");
    times[4] += PAPI_get_real_nsec() - t0;
    PAPI_stop(event_set, output_counters); 

    /************std::cout output for us that we don't want to profile***********/
    float FLOPS=0, TRAFFIC=0; 
    cout<<"\nFileReading:"<<endl;
    cout<<"Total elements read:"<<tot_nel_read<<", bytes:"<<tot_nel_read*sizeof(float)<<endl;
    output_profile(file_read_counters, times[0]);

    //profile processing the data
    cout<<"\nData reduction:"<<std::endl;
    output_profile(processing_counters, times[1]);

    cout<<"\nPk:"<<endl;
    output_profile(Pk_counters, times[2]);
      
    cout<<"\nPk serial:"<<endl;
    output_profile(Pk_counters_serial, times[3]);

    cout<<"\nOutput:"<<endl;
    output_profile(output_counters, times[4]);
    FLOPS = file_read_counters[0] + processing_counters[0] + Pk_counters[0] + output_counters[0];
    TRAFFIC = (file_read_counters[1] + processing_counters[1] +  Pk_counters[1] + output_counters[1])*sizeof(float);
    
    //calculate desired quantities    
    float Tp = static_cast<float>(times[0]+times[1]+times[2]+times[4]) * 1e-9;
    float Ts = static_cast<float>(times[0]+times[1]+times[3]+times[4]) * 1e-9;
    float Tp_sec = static_cast<float>(times[2]);
    float Ts_sec = static_cast<float>(times[3]);
    float rel_times[4];
    vector<int> it = {0,1,2,4};
    for(int i=0;i<4;++i)
        rel_times[i] = times[it[i]]/Tp*1e-7;
    size_t W = lst_Nr*lst_Nc*Nt*sizeof(float);  
    write_analysis(Pk_counters, Tp, rel_times, TRAFFIC, Tp_sec, Ts_sec, W, Ts, W, Tp,"results.txt");

    cout<<"\nTotals:\n";
    cout << "Execution time:       "<<Tp<<"s, File reading("<<rel_times[0]<<"%), Processing("<<rel_times[1]<<"%), Pk("<<rel_times[2]<<"%), Output("<<rel_times[3]<<"%)\n";     
    cout << "Traffic: (Bytes):     " << TRAFFIC << '\n';
    cout << "Flops:                " << FLOPS << '\n';
    cout << "Operational Intensity: " << FLOPS/tot_nel_read/sizeof(float) << '\n';
    
    //free(Pk);
    return 1;
}
