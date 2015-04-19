#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include <omp.h>

#ifndef BLOCK
#define BLOCK 16
#endif

extern "C" void  harris_opt(int  C, int  R, float * img, void * harris_void)
{

	float *dummy, *dummy1, *dummy2, *dummy3, *dummy4, *dummy5;
	float *Ixx, *Ixy, *Iyy, *Sxx, *Sxy, *Syy;

  float *harris = (float *)harris_void;
  
  #pragma omp parllel for
  for (int  ii = 0; ii < R; ii = (ii + BLOCK))
  {
		Ixx = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  Ixy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  Iyy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  Sxx = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
  	Sxy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  Syy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
    
		for (int  jj = 0; jj < C; jj = (jj + BLOCK))
		{
			int iblock = (R<ii+BLOCK)?(R-ii):BLOCK;
			int jblock = (C<jj+BLOCK)?(C-jj):BLOCK;
	  	
			for (int  i = 1; (i <= iblock); i = (i + 1))
			{
				for (int  j = 1; (j <= jblock); j = (j + 1))
				{
					float resx,resy;
					int index1 = (((-1 + (ii+i)) * (C + 2)) + (-1 + (jj+j)));		// [i-1] [j-1]
					int index2 = (((1 + (ii+i)) * (C + 2)) + (-1 + (jj+j))); 		// [i+1] [j-1]
					int index3 = (((-1 + (ii+i)) * (C + 2)) + (jj+j));			// [i-1] [j]
					int index4 = (((1 + (ii+i)) * (C + 2)) + (jj+j));				// [i+1] [j]
					int index5 = (((-1 + (ii+i)) * (C + 2)) + (1 + (jj+j))); 		// [i-1] [j+1]
					int index6 = (((1 + (ii+i)) * (C + 2)) + (1 + (jj+j)));		// [i+1] [j+1]
					int index7 = (((ii+i) * (C + 2)) + (-1 + (jj+j)));			// [i] [j-1]
					int index8 = (((ii+i) * (C + 2)) + (1 + (jj+j)));				// [i] [j+1]
					
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
										
					Ixx[((i * (2 + BLOCK)) + j)] = resx * resx;
					Iyy[((i * (2 + BLOCK)) + j)] = resy * resy;
					Ixy[((i * (2 + BLOCK)) + j)] = resx * resy;
 				}
			}
	  	
			for (int  i = 2; (i < iblock); i=i+1) {
				for (int  j = 2; (j < jblock); j=j+1) {
				
					int index1, index2, index3;
					index1 = (((-1 + i) * (2 + BLOCK)) + j);			// [i-1] [j]
					index2 = ((i * (2 + BLOCK)) + j);					// [i] [j]
					index3 = (((1 + i) * (2 + BLOCK)) + j);			// [i+1] [j]
	  
					Syy[index2] = Iyy[ index1-1 ] + 
                                 Iyy[ index1 ] + 
                                 Iyy[ index1+1 ] + 
                                 Iyy[ index2-1 ] + 
                                 Iyy[ index2 ] + 
                                 Iyy[ index2+1 ] + 
                                 Iyy[ index3-1 ] + 
                                 Iyy[ index3 ] + 
                                 Iyy[ index3+1 ];
				} 
			}
 
			for (int  i = 2; (i < iblock); i=i+1) {
				for (int  j = 2; (j < jblock); j=j+1){
				
					int index1, index2, index3;
					index1 = (((-1 + i) * (2 + BLOCK)) + j);			// [i-1] [j]
					index2 = ((i * (2 + BLOCK)) + j);					// [i] [j]
					index3 = (((1 + i) * (2 + BLOCK)) + j);			// [i+1] [j]
	  	
      		Sxy[index2] = Ixy[ index1-1 ] + 
                                 Ixy[ index1 ] +
                                 Ixy[ index1+1 ] + 
                                 Ixy[ index2-1 ] + 
                                 Ixy[ index2] + 
                                 Ixy[ index2+1 ] + 
                                 Ixy[ index3-1 ] + 
                                 Ixy[ index3 ] + 
                                 Ixy[ index3+1 ];
				} 
			}
			
			for (int  i = 2; (i < iblock); i=i+1) {
    		for (int  j = 2; (j < jblock); j=j+1) {
    		
				  int index1, index2, index3;
					index1 = (((-1 + i) * (2 + BLOCK)) + j);			// [i-1] [j]
					index2 = ((i * (2 + BLOCK)) + j);					// [i] [j]
					index3 = (((1 + i) * (2 + BLOCK)) + j);			// [i+1] [j]
	  
		      Sxx[index2] = Ixx[index1-1 ] + 
                                 Ixx[index1 ] + 
                                 Ixx[index1+1 ] + 
                                 Ixx[index2-1 ] + 
                                 Ixx[index2 ] + 
                                 Ixx[index2+1 ] + 
                                 Ixx[index3-1 ] + 
                                 Ixx[index3 ] + 
                                 Ixx[index3+1 ];
                                 
     
    		} 
			}
			
  		for (int  i = 2; (i < iblock); i++) {  
				for (int  j = 2; (j < jblock); j++) {
				
				  int index = (i * (2 + BLOCK)) + j;
				  int index2 = ((ii+i) * (2 + C)) + (jj+j);
				  float trace = 
          	Sxx[index] + Syy[index];

		      float det = 
	          Sxx[index] * Syy[index] - Sxy[index] * Sxy[index];

      		harris[index2] = det - (0.04f * trace * trace);
    		}
		  }
			
    }  
   	free(Ixx);
	  free(Iyy);
  	free(Ixy);
	  free(Sxx);
	  free(Sxy);
	  free(Syy);
  }
 }
