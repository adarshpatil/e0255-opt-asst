#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#include <algorithm>

#ifndef BLOCK
#define BLOCK 128
#endif

extern "C" void  harris_opt(int  C, int  R, float * img, void * harris_void)
{
	int iblock, jblock;
  float *harris = (float *)harris_void;
  
  #pragma omp parllel for
  for (int  ii = -1; ii <= R; ii = (ii + BLOCK))
  {
  	float Ixx [BLOCK+4] [BLOCK+4];
  	float Ixy [BLOCK+4] [BLOCK+4];
  	float Iyy [BLOCK+4] [BLOCK+4];
  	float Sxx [BLOCK+4] [BLOCK+4];
  	float Sxy [BLOCK+4] [BLOCK+4];
  	float Syy [BLOCK+4] [BLOCK+4];
    
		for (int  jj = -1; jj <= C; jj = (jj + BLOCK))
		{
			iblock = std::min(R,ii+BLOCK+3);
			jblock = std::min(C,jj+BLOCK+3);
			for (int  i = std::max(1,ii); i <= iblock; i = (i + 1))
			{
				#pragma GCC ivdep
				for (int  j = std::max(1,jj); j <= jblock; j = (j + 1))
				{
					float resx,resy;
					int index1 = (((-1 + i) * (C + 2)) + (-1 + j));		// [i-1] [j-1]
					int index2 = (((1 + i) * (C + 2)) + (-1 + j)); 		// [i+1] [j-1]
					int index3 = (((-1 + i) * (C + 2)) + j);			// [i-1] [j]
					int index4 = (((1 + i) * (C + 2)) + j);				// [i+1] [j]
					int index5 = (((-1 + i) * (C + 2)) + (1 + j)); 		// [i-1] [j+1]
					int index6 = (((1 + i) * (C + 2)) + (1 + j));		// [i+1] [j+1]
					int index7 = ((i * (C + 2)) + (-1 + j));			// [i] [j-1]
					int index8 = ((i * (C + 2)) + (1 + j));				// [i] [j+1]
					
					// X derivative
					resx = (img[index1] * -0.0833333333333f) + 
										(img[index2] * 0.0833333333333f) + 
										(img[index3] * -0.166666666667f) + 
										(img[index4] * 0.166666666667f) + 
										(img[index5] * -0.0833333333333f) + 
										(img[index6] * 0.0833333333333f);
                                
					// Y derivative
					resy = (img[index1] * -0.0833333333333f) + 
										(img[index5] * 0.0833333333333f) + 
										(img[index7] * -0.166666666667f) + 
										(img[index8] * 0.166666666667f) + 
										(img[index2] * -0.0833333333333f) + 
										(img[index6] * 0.0833333333333f);
										
					Ixx[i-ii][j-jj] = resx * resx;
					Iyy[i-ii][j-jj] = resy * resy;
					Ixy[i-ii][j-jj] = resx * resy;
 				}
			}
	  	
	  	iblock = std::min(R-1,ii+BLOCK+2);
	  	jblock = std::min(C-1,jj+BLOCK+2);
	  	
			for (int  i = std::max(2,ii+1); i <= iblock; i=i+2) {
				#pragma GCC ivdep
				for (int  j = std::max(2,jj+1); j <= jblock; j=j+1) {
	  
					Syy[i-ii][j-jj] = Iyy[i-ii-1][j-jj-1] + 
                                 Iyy[i-ii-1][j-jj] + 
                                 Iyy[i-ii-1][j-jj+1] + 
                                 Iyy[i-ii][j-jj-1] + 
                                 Iyy[i-ii][j-jj] + 
                                 Iyy[i-ii][j-jj+1] + 
                                 Iyy[i-ii+1][j-jj-1] + 
                                 Iyy[i-ii+1][j-jj] + 
                                 Iyy[i-ii+1][j-jj+1];
                                 
					Syy[i-ii+1][j-jj] = Iyy[i-ii][j-jj-1] + 
                                 Iyy[i-ii][j-jj] + 
                                 Iyy[i-ii][j-jj+1] + 
                                 Iyy[i-ii+1][j-jj-1] + 
                                 Iyy[i-ii+1][j-jj] + 
                                 Iyy[i-ii+1][j-jj+1] + 
                                 Iyy[i-ii+2][j-jj-1] + 
                                 Iyy[i-ii+2][j-jj] + 
                                 Iyy[i-ii+2][j-jj+1];                                 
				} 
			}
 
			for (int  i = std::max(2,ii+1); i <= iblock; i=i+2) {
				#pragma GCC ivdep
				for (int  j = std::max(2,jj+1); j <= jblock; j=j+1){
				
      		Sxy[i-ii][j-jj] = Ixy[i-ii-1][j-jj-1] + 
                                 Ixy[i-ii-1][j-jj] + 
                                 Ixy[i-ii-1][j-jj+1] + 
                                 Ixy[i-ii][j-jj-1] + 
                                 Ixy[i-ii][j-jj] + 
                                 Ixy[i-ii][j-jj+1] + 
                                 Ixy[i-ii+1][j-jj-1] + 
                                 Ixy[i-ii+1][j-jj] + 
                                 Ixy[i-ii+1][j-jj+1];
                                 
      		Sxy[i-ii+1][j-jj] = Ixy[i-ii][j-jj-1] + 
                                 Ixy[i-ii][j-jj] + 
                                 Ixy[i-ii][j-jj+1] + 
                                 Ixy[i-ii+1][j-jj-1] + 
                                 Ixy[i-ii+1][j-jj] + 
                                 Ixy[i-ii+1][j-jj+1] + 
                                 Ixy[i-ii+2][j-jj-1] + 
                                 Ixy[i-ii+2][j-jj] + 
                                 Ixy[i-ii+2][j-jj+1];                                 
				} 
			}
			
			for (int  i = std::max(2,ii+1); (i <= iblock); i=i+2) {
				#pragma GCC ivdep
    		for (int  j = std::max(2,jj+1); (j <= jblock); j=j+1) {
    			
		      Sxx[i-ii][j-jj] = Ixx[i-ii-1][j-jj-1] + 
                                 Ixx[i-ii-1][j-jj] + 
                                 Ixx[i-ii-1][j-jj+1] + 
                                 Ixx[i-ii][j-jj-1] + 
                                 Ixx[i-ii][j-jj] + 
                                 Ixx[i-ii][j-jj+1] + 
                                 Ixx[i-ii+1][j-jj-1] + 
                                 Ixx[i-ii+1][j-jj] + 
                                 Ixx[i-ii+1][j-jj+1];
                                 
		      Sxx[i-ii+1][j-jj] = Ixx[i-ii][j-jj-1] + 
                                 Ixx[i-ii][j-jj] + 
                                 Ixx[i-ii][j-jj+1] + 
                                 Ixx[i-ii+1][j-jj-1] + 
                                 Ixx[i-ii+1][j-jj] + 
                                 Ixx[i-ii+1][j-jj+1] + 
                                 Ixx[i-ii+2][j-jj-1] + 
                                 Ixx[i-ii+2][j-jj] + 
                                 Ixx[i-ii+2][j-jj+1];                                 
     
    		} 
			}
			
			iblock = std::min(R-1,ii+BLOCK+1);
			jblock = std::min(C-1,jj+BLOCK+1);
			
  		for (int  i = ii+2; (i <= iblock); i++) {
  			#pragma GCC ivdep 
				for (int  j = jj+2; (j <= jblock); j++) {
				
				  int index = (i * (2 + C)) + j;
				  float trace = 
          	Sxx[i-ii][j-jj] + Syy[i-ii][j-jj];

		      float det = 
	          Sxx[i-ii][j-jj] * Syy[i-ii][j-jj] - Sxy[i-ii][j-jj] * Sxy[i-ii][j-jj];

      		harris[index] = det - (0.04f * trace * trace);
    		}
		  }
    }
  }
 }
