#include "Geotiff.hpp"
#include <algorithm>
#include "papi.h"
#include <numeric>
#include <filesystem>
#include <regex>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>

//Just a wrapper to index 2D image stored as flat array
class Image
{
    private:
        int Nrow;
        int Ncol;
       
    public:
        vector<float> vec;
        float *data;
        
        Image(vector<float> &_data, int _dims[2])
        {
            //get image dimensions
            Nrow = _dims[0];
            Ncol = _dims[1];
            data = _data.data();    
            vec = _data;
        }

        float &operator()(int row, int col)
        {
            return data[row+col*Nrow]; //col major order (to match eigen)
        }
};

void write_to_file(map<int, map<int,float>> maxes, string filename)
{
    ofstream ofile(filename.c_str()); 
    vector<string> months= {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    size_t ncols;
    for(auto y : maxes){
        ncols = y.second.size();   
        break;
    }
    for(size_t i=0;i<ncols;i++)
        ofile << months[i] << " ";
    ofile<<endl; 
    for( auto y : maxes){
        ofile << y.first << " ";
        for(auto m : y.second)
            ofile << m.second <<" ";
        ofile<<endl;
    }
    ofile.close();
}

void write_binary(float* data, int N, int Npc, string fname)
{
    ofstream file(fname.c_str(), ios::out|ios::binary);
    file.write((char*)&N, sizeof(int));
    file.write((char*)&Npc, sizeof(int));
    file.write((char*)data, N*Npc*sizeof(float));
    file.close();
}


//just organizes images by date
struct  DateContainer
{
   map<int, map<int, map<int, vector<float>>>> images;
   size_t size = 0;

   //accesses image(s) for particular year/month/day
   vector<float> &operator()(int year, int month, int day){
        return images[year][month][day];
   }

   //accesses image(s) for particular year/month/day
   vector<float> &operator()(string year, string month, string day){
        return images[atoi(year.c_str())][atoi(month.c_str())-1][atoi(day.c_str())];
   }

   //add image while keeping count of how many you've added
   void add(string year, string month, string day, vector<float> & addv)
   {
        images[atoi(year.c_str())][atoi(month.c_str())-1][atoi(day.c_str())] = addv;
        size++;
        return;    
   }   
};

void output_profile( long long int counters[4], long long int dt)
{
    long long sflops = counters[0];      
    long long dflops = counters[1]; 
    long long traffic = counters[2]*sizeof(float); 
    long long total_l1d_cache_misses = counters[3];  // number of such instructions
    
    double twall = static_cast<double>(dt) * 1.0e-9; // seconds
    cout << "Wall time(s):                 " << twall << '\n';
    cout << "Flops (Single prec):          " << sflops << '\n';
    cout << "Flops (Double prec):          " << dflops << '\n';
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
    int event_set[2] = {PAPI_NULL, PAPI_NULL};
    int events[4]    = {PAPI_SP_OPS, PAPI_DP_OPS, PAPI_LST_INS, PAPI_L1_DCM};
    //int memevents[1] = {PAPI_LST_INS};
    long long int file_read_counters[4]  = {0,0,0,0};
    long long int processing_counters[4] = {0,0,0,0};
    long long int ssa_counters[4] = {0,0,0,0};
    long long int output_counters[4] = {0,0,0,0};
    //long long int total_lst_ins=0;
 
    PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&event_set[0]);
    //PAPI_create_eventset(&event_set[1]);
    PAPI_add_events(event_set[0], events, 4);
    //PAPI_add_events(event_set[1], memevents, 4); 
    long long int t0;
    long long int times[4] = {0,0,0,0};

    //profile reading the files
    map<int, map<int,float>> maxes;
    map<int, map<int,float>> mins;
    map<int, map<int,float>> means;

    /********************************read the files*************************/
    DateContainer data;
    size_t tot_nel_read=0;
    t0 = PAPI_get_real_nsec();
    PAPI_start(event_set[0]);
    for(const auto& pszFname : filesystem::recursive_directory_iterator(pszDir))
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
        data.add(year, month, day, gtiff.GetRasterBand(1));
        tot_nel_read += data(year, month, day).size();
        //cout<<"tot_nel_read:"<<tot_nel_read<<endl;
    }   
    //PAPI_stop(event_set[0], file_read_counters);
    PAPI_accum(event_set[0], file_read_counters);
    times[0] += PAPI_get_real_nsec() - t0;
    
    //cout<<"starting processing: N:"<<data.size<<endl;
    /***************************do processing******************************/
    //first need to convert counts to actual data values and throw out boundary values
    //form time series for total NDVI while we do this.
    vector<float> tot_ndvi(data.size); 
    //loop over images 
    vector<tuple<int,int,int>> dates;
    size_t cnt=0;
    t0 = PAPI_get_real_nsec();
    //PAPI_start(event_set[0]);
    for(auto &yearv : data.images){
        int year = yearv.first;
        for(auto &monthv : yearv.second){
            int month = monthv.first;
            float max = -numeric_limits<float>::max();
            float min =  numeric_limits<float>::min();
            float mean = 0, dc=0;
            for(auto &dayv : monthv.second){
                int day = dayv.first;
                float _totndvi = 0, pc=0;
                for(float &val : dayv.second){
                    if(val<255){
                        float tmp = (val-125)/125.0;
                        _totndvi += tmp;
                        max = tmp>max ? tmp : max;
                        min = tmp<min ? tmp : min;
                        mean += tmp; 
                        dc++; pc++;               
                    }
                }
                tot_ndvi[cnt] = _totndvi/pc;
                dates.push_back(make_tuple(year,month,day));
                cnt++;
            }
            maxes[year][month] = max > maxes[year][month] ? max : maxes[year][month];
            mins[year][month]  = min < mins[year][month] ? min : mins[year][month];
            means[year][month] = mean/dc;  
        }
    }
    PAPI_accum(event_set[0], processing_counters);
    times[1] += PAPI_get_real_nsec() - t0;
    //if(cnt!=tot_ndvi.size()) throw(runtime_error("map loop doesn't match files read!\n")); //check we don't want to profile    
 
    /*********************************now do the SSA*****************************/
    t0 = PAPI_get_real_nsec();
    //PAPI_start(event_set[0]);
    //form trajectory matrix
    int N = tot_ndvi.size();
    int L = int(N/2);
    int K = N-L+1;
    Eigen::MatrixXf mat(L, K);
    for(int col=0;col<K;col++){
        for(int row=0;row<L;row++){
            mat(row, col) = tot_ndvi[row+col];
        }
    }
    //now do svd
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(mat,Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXf U = svd.matrixU();
    Eigen::MatrixXf V = svd.matrixV();
    Eigen::VectorXf s = svd.singularValues();

    // now generate resulting time series
    int LKmin = min(L, N-L+1), LKmax = max(L, N-L+1);
    int Npc = s.size();
    vector<float> results(N*Npc);
    for(int pc=0; pc<Npc; pc++){
        //reconstruct pc from U(:,pc)*s(pc)*V(:,pc)^T
        float sv = s(pc);
        for(int col=0;col<K;col++){
            float v = V(pc,col);
            for(int row=0; row<L; row++)
                results[pc*N+row+col] = sv*U(row,pc)*v;
        }
        //do anti-diagonal averaging
        for(int l=0;l<LKmin;l++)
            results[pc*N+l] /= l+1;
        for(int l=LKmin; l<LKmax; l++)
            results[pc*N+l] /= LKmin;
        for(int l=LKmax; l<N;l++)
            results[pc*N+l] /= (N-l);
    }
    times[2] += PAPI_get_real_nsec() - t0;
    PAPI_accum(event_set[0], ssa_counters);

    /*********************************now do the output**************************/
    t0 = PAPI_get_real_nsec();

    //Write out a text file which has timing data so python can use it
    write_to_file(maxes, "maxes.txt");
    write_to_file(means, "means.txt");
    write_to_file(mins, "mins.txt");  
    write_binary(results.data(), N, Npc, "ssa");
    times[3] += PAPI_get_real_nsec() - t0;
    PAPI_stop(event_set[0], output_counters); 

    /************std::cout output for us that we don't want to profile***********/
    /*cout<<"\nNDVI trends:"<<endl;
    for(auto &M : maxes){
        int year = M.first;
        cout<<"Year:"<<year<<"\n";
        for(auto &m : M.second){
            cout<<"Month:"<<m.first<<", max ndvi: "<<m.second<<", mean:"<<means[year][m.first]<<", min:"<<mins[year][m.first]<<endl;    
        }
    }*/

    double FLOPS=0, TRAFFIC=0; 
    std::cout<<"\nFileReading:"<<std::endl;
    cout<<"Total elements read:"<<tot_nel_read<<", bytes:"<<tot_nel_read*sizeof(float)<<endl;
    output_profile(file_read_counters, times[0]);

    //profile processing the data
    std::cout<<"\nData reduction:"<<std::endl;
    output_profile(processing_counters, times[1]);

    cout<<"\nSSA:"<<endl;
    output_profile(ssa_counters, times[2]);
      
    cout<<"\nOutput:"<<endl;
    output_profile(output_counters, times[3]);
    
    FLOPS = file_read_counters[0] + processing_counters[0] + ssa_counters[0] + output_counters[0];
    TRAFFIC = (file_read_counters[2] + output_counters[2])*sizeof(float);
    double total_time = static_cast<double>(times[0]+times[1]+times[2]+times[3]) * 1e-9;
    double rel_times[4];
    for(int i=0;i<4;i++)
        rel_times[i] = times[i]/total_time*1e-7;

    cout<<"\nTotals:\n";
    cout << "Execution time:"<<total_time<<"s, File reading("<<rel_times[0]<<"%), Processing("<<rel_times[1]<<"%), SSA("<<rel_times[2]<<"%), Output("<<rel_times[3]<<"%)\n";     
    //cout << "Traffic1: (Bytes):             " << total_lst_ins*sizeof(float) << '\n';
    cout << "Traffic: (Bytes):             " << TRAFFIC << '\n';
    cout << "Flops:                        " << FLOPS << '\n';
    cout << "Opertional Intensity:         " << FLOPS/tot_nel_read/sizeof(float) << '\n';
    return 1;
}
