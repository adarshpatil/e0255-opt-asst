#include <stdio.h>
#include <stdlib.h>

// check L1 / L2 size and tune
#ifndef BLOCK
#define BLOCK 128
#endif

//check 4k aliasing and tune
#ifndef PAD
#define PAD 1024
#endif

void  harris_opt(int  C, int  R, float * img, float *& harris)
{

	float *dummy, *dummy1, *dummy2, *dummy3, *dummy4, *dummy5;
	float *Ixx, *Ixy, *Iyy, *Sxx, *Sxy, *Syy;

  harris = (float *) (malloc((sizeof(float ) * ((2 + R) * (2 + C)))));	
  //dummy = (float *) (malloc(PAD));
  #pragma omp parllel for
  for (int  ii = 0; ii < R; ii = (ii + BLOCK))
  {
		Ixx = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  //dummy1 = (float *) (malloc(PAD));
	  Ixy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  //dummy2 = (float *) (malloc(PAD));
	  Iyy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  //dummy3 = (float *) (malloc(PAD));
	  Sxx = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
  	//dummy4 = (float *) (malloc(PAD));
  	Sxy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
	  //dummy5 = (float *) (malloc(PAD));
	  Syy = (float *) (malloc((sizeof(float ) * ((2 + BLOCK) * (2 + BLOCK)))));
    
		for (int  jj = 0; jj < C; jj = (jj + BLOCK))
		{
			int iblock = std::min(R,ii+BLOCK) - ii;
			int jblock = std::min(C,jj+BLOCK) - jj;
	  	
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
	  //free(dummy); 
  	//free(dummy1); free(dummy2); free(dummy3); free(dummy4); free(dummy5);
  }
 }
