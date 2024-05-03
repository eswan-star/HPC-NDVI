#include "transforms.hpp"

#define PI2 6.28318530717958

#define set_mm512(cos,l,n,N)\
    _mm512_set_ps(cos((l)*(n+15)*PI2/N),cos((l)*(n+14)*PI2/N),cos((l)*(n+13)*PI2/N),\
		  cos((l)*(n+12)*PI2/N),cos((l)*(n+11)*PI2/N),cos((l)*(n+10)*PI2/N),\
                  cos((l)*(n+9)*PI2/N),cos((l)*(n+8)*PI2/N),cos((l)*(n+7)*PI2/N),\
                  cos((l)*(n+6)*PI2/N),cos((l)*(n+5)*PI2/N),cos((l)*(n+4)*PI2/N),\
                  cos((l)*(n+3)*PI2/N),cos((l)*(n+2)*PI2/N),cos((l)*(n+1)*PI2/N),\
	          cos((l)*(n)*PI2/N));

#define PIPL_CHUNK 7
#define SIMD_CHUNK 16

#define GET_F_COL(cos,l,n,N)\
    f_col[0] = set_mm512(cos,l,n,N)\
    f_col[1] = set_mm512(cos,l+1,n,N)\
    f_col[2] = set_mm512(cos,l+2,n,N)\
    f_col[3] = set_mm512(cos,l+3,n,N)\
    f_col[4] = set_mm512(cos,l+4,n,N)\
    f_col[5] = set_mm512(cos,l+5,n,N)\
    f_col[6] = set_mm512(cos,l+6,n,N)

#define ACCUMULATE()\
    accum[0] = _mm512_fmadd_ps(img_slice, f_col[0], accum[0]);\
    accum[1] = _mm512_fmadd_ps(img_slice, f_col[1], accum[1]);\
    accum[2] = _mm512_fmadd_ps(img_slice, f_col[2], accum[2]);\
    accum[3] = _mm512_fmadd_ps(img_slice, f_col[3], accum[3]);\
    accum[4] = _mm512_fmadd_ps(img_slice, f_col[4], accum[4]);\
    accum[5] = _mm512_fmadd_ps(img_slice, f_col[5], accum[5]);\
    accum[6] = _mm512_fmadd_ps(img_slice, f_col[6], accum[6]);

#define F(op,m,k,N)\
	op(PI2*m*k/N)

using namespace std;

void right_transform(float * __restrict__ img_data, const int &Nt, const int &Nr, const int &Nc)
{
 
    alignas(64) float row_tmp[Nc];
    alignas(64) float img_data_tmp[Nc];
    __m512 img_slice, mask;
    __m512 f_col[PIPL_CHUNK], accum[PIPL_CHUNK]; 
    int lstop = floor((float)Nc/PIPL_CHUNK)*PIPL_CHUNK;
    int nstop = floor((float)Nc/SIMD_CHUNK)*SIMD_CHUNK; 
    for(int T=0; T<Nt; ++T)
    {
	    for(int m=0; m<Nr; m++)
	    { 
		    //zero temporary memory used for storing row for result 
		    memset(row_tmp, 0, Nc*sizeof(float));
		    memcpy(img_data_tmp, &img_data[T*Nc*Nr + m*Nc], Nc*sizeof(float));
		    for(int l=0; l<lstop; l+=PIPL_CHUNK)
		    {  
			//zero accumulators
			for(int i=0; i<PIPL_CHUNK; ++i)
			    accum[i] = _mm512_set1_ps(0.0);
			//do chunked loop
			for(int n=0; n<nstop; n+=SIMD_CHUNK)
			{
			    //load image slice
			    img_slice = _mm512_load_ps(&img_data_tmp[n]);
			    
			    //get chunk of transform-matrix
			    GET_F_COL(cos,l,n,Nc)
				    
			    //multiply with row and store
			    ACCUMULATE()
			}
			
			//Do remainder
			int rem = Nc - nstop;
			alignas(64) float msk[SIMD_CHUNK];
			memset(msk,0,SIMD_CHUNK*sizeof(float));
			for(int m=0;m<rem;m++)
			    msk[m] = 1;
			
			//This will read past the end of the matrix we're currently
			//considering, but those entries wil be zeroed out. Would be
			//a segfault though if was last row of matrix (fortunately, 
			//never will be because that's handled outside this main loop.)
			mask = _mm512_load_ps(msk); 
			img_slice = _mm512_load_ps(&img_data_tmp[nstop]);
			img_slice = _mm512_mul_ps(img_slice, mask);
			GET_F_COL(cos,l,nstop,Nc)
			ACCUMULATE()	
			//do final horizontal add and write to dest 
			alignas(64) float tmp[SIMD_CHUNK];
			for(int li=0;li<PIPL_CHUNK;++li){
			    _mm512_store_ps(tmp, accum[li]);
			    row_tmp[l+li] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+
				tmp[6]+tmp[7]+tmp[8]+tmp[9]+tmp[10]+tmp[11]+tmp[12]+
				tmp[13]+tmp[14]+tmp[15];
			}	
		    }
		    //now need to do all of the above but when there's less than PIPL_CHUNK columns left
		    //for now do in loop b/c can't know remaining dimensions until runtime.  
		    for(int l=lstop;l<Nc;l++){
			accum[0] = _mm512_set1_ps(0.0);
			for(int n=0;n<nstop;n+=SIMD_CHUNK){
			    img_slice = _mm512_load_ps(&img_data_tmp[n]);
			    f_col[0] = set_mm512(cos,l,n,Nc)
			    accum[0] = _mm512_fmadd_ps(img_slice, f_col[0], accum[0]);
			}
			//Now can't do vectorized read because would go past end of matrix;
			//so loop:
			alignas(64) float tmp[SIMD_CHUNK];
			_mm512_store_ps(tmp, accum[0]);
			row_tmp[l] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+
				tmp[6]+tmp[7]+tmp[8]+tmp[9]+tmp[10]+tmp[11]+tmp[12]+
				tmp[13]+tmp[14]+tmp[15]; 
			for(int n=nstop;n<Nc;n++){
			    //now need to ad up last un-accumulated fourier terms
			    row_tmp[l] += F(cos,l,n,Nc) * img_data_tmp[n];
			}
		    }
		    //congrats, you've finished a row. Now write it out.
		    memcpy(&img_data[T*Nr*Nc + m*Nc], row_tmp, Nc*sizeof(float));
	    }
    } 
    return;
}


void left_transform_serial(float* __restrict__ img_data, const int &Nt, const int&Nr, const int &Nc)
{
    float mat_tmp[Nr*Nc];// = (float*)malloc(Nr*Nc*sizeof(float));
    memset(mat_tmp,0,sizeof(float)*Nc*Nr);
    int chunk = 32768/(2.0*sizeof(float)); //ACJ chunk size to fit in cache
    //ACJ: outer loop for time
    for(int T=0; T<Nt; T++)
    {
        //#pragma omp parallel reduction(max:dt)
        { 
            //#pragma omp for 
            for(size_t m=0;m<Nr;m++)
            {   
                int start = 0;
                int stop = min(start+chunk, Nc);
                while(start<Nc)
                {
                    for(size_t k=0;k<Nr;k++)
                    {
                        float f_mk = F(cos,m,k,Nr);
                        for(size_t n=start;n<stop;n++)
                        {
                            mat_tmp[m*Nc+n] += f_mk*img_data[T*Nr*Nc + m*Nc + n];
                        } 
                    }
                start = stop;
                stop = min(start+chunk, Nc);
                }
            }
        }
	memcpy(&img_data[T*Nr*Nc], mat_tmp, Nr*Nc*sizeof(float));
	memset(mat_tmp,0,Nc*Nr*sizeof(float));
    }
}

//ACJ: for now just gets even term but will add even/odd combinations after being completely 
//done profiling.
void get_Pk_kernel(float*__restrict__ img_data, const int &Nt, const int &Nr, const int &Nc)
{
   //This is pretty simple for now, just apply one transform then the other. 
   cout<<"Right transform."<<endl;
   right_transform(img_data, Nt, Nr, Nc);
   cout<<"Left transform."<<endl;
   left_transform_serial(img_data, Nt, Nr, Nc);
}

void get_Pk_kernel_serial(float*__restrict__ img_data, const int &Nt, const int &Nr, const int &Nc)
{
   right_transform(img_data, Nt, Nr, Nc);
   left_transform_serial(img_data, Nt, Nr, Nc); 
}

