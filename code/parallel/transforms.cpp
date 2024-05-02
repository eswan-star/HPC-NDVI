#include <algorithm>
#include "papi.h"
#include <numeric>
#include <fstream>
#include <x86intrin.h>
#include <stdalign.h>
#include <cstring>


#define PI2 6.28318530717958

#define set_mm512(cos,l,n)\
    _mm512_set_ps(cos((l)*(n+15)*PI2/Nc),cos((l)*(n+14)*PI2/Nc),cos((l)*(n+13)*PI2/Nc),\
		  cos((l)*(n+12)*PI2/Nc),cos((l)*(n+11)*PI2/Nc),cos((l)*(n+10)*PI2/Nc),\
                  cos((l)*(n+9)*PI2/Nc),cos((l)*(n+8)*PI2/Nc),cos((l)*(n+7)*PI2/Nc),\
                  cos((l)*(n+6)*PI2/Nc),cos((l)*(n+5)*PI2/Nc),cos((l)*(n+4)*PI2/Nc),\
                  cos((l)*(n+3)*PI2/Nc),cos((l)*(n+2)*PI2/Nc),cos((l)*(n+1)*PI2/Nc),\
	          cos((l)*(n)*PI2/Nc));

#define PIPL_CHUNK 7
#define SIMD_CUNK 16

#define GET_F_COL\
    f_col[0] = set_mm512(cos,l,n)\
    f_col[1] = set_mm512(cos,l+1,n)\
    f_col[2] = set_mm512(cos,l+2,n)\
    f_col[3] = set_mm512(cos,l+3,n)\
    f_col[4] = set_mm512(cos,l+4,n)\
    f_col[5] = set_mm512(cos,l+5,n)\
    f_col[6] = set_mm512(cos,l+6,n)

#define ACCUMULATE\
    accum[0] = _mm512_fmadd(img_slice, f_col[0], accum[0]);\
    accum[1] = _mm512_fmadd(img_slice, f_col[1], accum[1]);\
    accum[2] = _mm512_fmadd(img_slice, f_col[2], accum[2]);\
    accum[3] = _mm512_fmadd(img_slice, f_col[3], accum[3]);\
    accum[4] = _mm512_fmadd(img_slice, f_col[4], accum[4]);\
    accum[5] = _mm512_fmadd(img_slice, f_col[5], accum[5]);\
    accum[6] = _mm512_fmadd(img_slice, f_col[6], accum[6]);

void inline right_transform(float *img_data, int &Nr, int &Nc, int &Nt)
{
    float row_tmp[Nc];
    __mm512 img_slice, f_col[PIPL_CHUNK], accum[PIPL_CHUNK], mask; 
    int lstop = floor(Nc/PIPL_CHUNK);
    int nstop = floor(Nc/SIMD_CHUNK);
    for(int m=0; m<Nr; m++)
    { 
	    //zero temporary memory used for storing row fo result 
	    memset(row_tmp,0,sizeof(float)*Nc);
	    for(l=0;l<lstop;l+=PIPL_CHUNK)
	    {   
		//zero accumulators
		for(int i=0;i<PIPL_CHUNK;++i)
			accum[i] = _mm512_set1_ps(0.0);
		//do chunked loop
		for(int n=0;n<nstop; n+=SIMD_CHUNK)
		{
		    //load image slice
		    img_slice = _mm512_load_ps(img_data[m*Nc+n]);
		    //get chunk of transform-matrix
		    GET_F_COL
			    
		    //multiply with row and store
		    ACCUMULTATE
		}

		//do remainder
		float msk[SIMD_CHUNK];
		memset(msk,0,SIMD_CHUNK*sizeof(float));
		for(int m=0;m<=rem;m++)
		    msk[m] = 1;
		//this will read past the end of the matrix we're currently
		//considering, but those entries wil be zeroed out. Would be
		//a segfault though if was last row of matrix (fortunately, 
		//never will be because that's handled outside this main loop.)
		mask = _mm512_load_ps(msk); 
		img_slice =_mm512_load_ps(img_data[m*Nc+n]);
		img_slice = _mm512_mult_ps(img_slice, mask);
		GET_F_COL
		ACCUMULATE
		
		//do final horizontal add and write to dest 
		float tmp[SIMD_CHUNK];
		for(int li=0;li<PIPL_CHUNK;++li){
		    _mm512_store_ps(tmp, accum[li]);
		    row_tmp[l+li] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+
			tmp[6]+tmp[7]+tmp[8]+tmp[9]+tmp[10]+tmp[11]+tmp[12]+
			tmp[13]+tmp[14]+tmp[15];
		}
		
	    }
	    //now need to do all of the above but when there's less than PIPL_CHUNK columns left
	    //for now do in loop b/c can't know remaining dimensions until runtime. 
	    for(l=lstop;l<Nc;l++){
		accum[0] = _mm512_set1_ps(0.0);
		for(int n=0;n<nstop;n+=SIMD_CHUNK){
		    img_slice = _mm512_load_ps(img_data[m*Nc+n]);
		    f_col[0] = set_mm512(cos,l,n)
		    accum[0] = _mm512_fmadd(img_slice, f_col[0], accum[0]);
		}
		//Now can't do vectorized read because would go past end of matrix;
		//so loop:
		float f[SIMD_CHUNK];
		float tmp[SIMD_CHUNK];
		f_col[0] = set_mm512(cos,l,n)
	        _mm512_store_ps(f, f_col[0]);//just to get fourier terms
	        _mm512_store_ps(tmp, accum);
	        row_tmp[l] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+
			tmp[6]+tmp[7]+tmp[8]+tmp[9]+tmp[10]+tmp[11]+tmp[12]+
			tmp[13]+tmp[14]+tmp[15]; 
	        for(n=nstop;n<Nc;n++){
		    //now need to ad up last un-accumulated fourier terms
		    row_tmp[l] += f[n]*img_data[m*Nc+n];
		}
	    }
	    //congrats, you've finished a row. Now write it out.
 	    memcpy(img_data[m*Nc], row_tmp, Nc*sizeof(float));
    } 
    return;
}


void inline left_transform()
{

}

