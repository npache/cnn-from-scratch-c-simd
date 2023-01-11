//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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include "../tp2.h"
#include "tiempo.h"
#include "libbmp.c"
#include "utils.c"
#include "imagenes.c"
#include "../cli.c"
#include "multiplicarMatrices.h"
int epch_counter=-1;

int batch_counter=-1;

void matrix_div(Matrix_t matrix, Matrix_t matrix_B){
  for(int i=0;i<matrix.height;i++){
    for(int j=0;j<matrix.width;j++){
      matrix.matrix[i][j]=matrix.matrix[i][j]/matrix_B.matrix[i][j];
    }
  }
}
void print_filter(char* s, Filter_t filt){
FILE * fp;
fp = fopen (s, "w+");
for(int i=0;i<filt.depth_out;i++){
      for(int j=0;j<filt.depth_in;j++){
        for(int k=0;k<filt.f;k++){
          for(int l=0;l<filt.f;l++){
    fprintf(fp,"%.10Lf \n",filt.filter[i][j].matrix[k][l]);
    }
  } 
}
}
fclose(fp);
}

void print_image(char* s, Image_t m){
FILE * fp;
fp = fopen (s, "w+");
for(int i = 0;i< m.channels;i++){
for(int j = 0;j< m.height;j++){
  for(int k = 0;k< m.width;k++){
    fprintf(fp,"%.10Lf \n",m.image[i].matrix[j][k]);
    }
  } 
}
fclose(fp);
}


void print_matrix(char* s, Matrix_t m){
FILE * fp;
fp = fopen (s, "w+");
for(int i = 0;i< m.height;i++){
  for(int j = 0;j< m.width;j++){
    fprintf(fp,"%.10Lf \n",m.matrix[i][j]);
    }
  } 
fclose(fp);
}

void copy_filter(Filter_t tot_f1_copy,Filter_t tot_f1){
  for(int i=0;i<tot_f1_copy.depth_out;i++){
    for(int j=0;j<tot_f1_copy.depth_in;j++){
      for(int k=0;k<tot_f1_copy.f;k++){
        for(int l=0;l<tot_f1_copy.f;l++){
        tot_f1_copy.filter[i][j].matrix[k][l]=tot_f1.filter[i][j].matrix[k][l];

        } 
      }
    }
  }

}

void copy_matrix(Matrix_t tot_f1_copy,Matrix_t tot_f1){
  for(int i=0;i<tot_f1_copy.height;i++){
    for(int j=0;j<tot_f1_copy.width;j++){
 
        tot_f1_copy.matrix[i][j]=tot_f1.matrix[i][j];

        } 
      }
    }

void matrix_initialize_n(Matrix_t* matrix, uint32_t height, uint32_t width){
  matrix->height=height;
  matrix->width=width;
  matrix->matrix = (long double **)malloc(height*sizeof(long double*));
  
  for(int i = 0;i< height;i++){
      
    matrix->matrix[i]= (long double *)malloc( width*sizeof(long double));
    for(int j = 0;j< width;j++){
      matrix->matrix[i][j]=0;
      
    }
  }
      
}

void matrix_initialize_n_free(Matrix_t matrix){
  
  for(int i = 0;i< matrix.height;i++){
    free(matrix.matrix[i]);
  }
   free(matrix.matrix);
}
void matrix_mul_scalar(Matrix_t matrix, long double x){
  for(int i=0;i<matrix.height;i++){
    for(int j=0;j<matrix.width;j++){
      matrix.matrix[i][j]=matrix.matrix[i][j]*x;
    }
  }
}

void matrix_add_scalar(Matrix_t matrix, long double x){
  for(int i=0;i<matrix.height;i++){
    for(int j=0;j<matrix.width;j++){
      matrix.matrix[i][j]=matrix.matrix[i][j]+x;
    }
  }
}

void matrix_sub_scalar(Matrix_t matrix, long double x){
  for(int i=0;i<matrix.height;i++){
    for(int j=0;j<matrix.width;j++){
      matrix.matrix[i][j]=matrix.matrix[i][j]-x;
    }
  }
}
void matrix_div_scalar(Matrix_t matrix, long double x){
  for(int i=0;i<matrix.height;i++){
    for(int j=0;j<matrix.width;j++){
      matrix.matrix[i][j]=matrix.matrix[i][j]/x;
    }
  }
}

void matrix_sqrt(Matrix_t matrix,int debugger){
  for(int i=0;i<matrix.height;i++){
    for(int j=0;j<matrix.width;j++){
      matrix.matrix[i][j]=sqrt(matrix.matrix[i][j]);
    }
  }
}
     


void dot_n(Matrix_t matrix_A,Matrix_t matrix_B,Matrix_t res)
{
    assert(matrix_A.width==matrix_B.height && "The dimension of the matrices do not match");

    for(int i = 0;i< matrix_A.height;i++){
        for(int j = 0;j< matrix_B.width;j++){
            long double sum = 0;
            for(int k = 0;k< matrix_A.width;k++){

                sum+= matrix_A.matrix[i][k]*matrix_B.matrix[k][j];
            }   
            res.matrix[i][j]=sum;  
            
                
        }
    }

}



void matrices_sub(
    long double **matrix_A,
    long double **matrix_B,
    int width_A,
    int height_A,
    int width_B,
    int height_B){
    assert(width_A==width_B && "The dimension of the matrices do not match");
    assert(height_A==height_B && "The dimension of the matrices do not match");
    for(int i = 0;i< height_A;i++){
        for(int j = 0;j< width_A;j++){
          matrix_A[i][j] = matrix_A[i][j]-matrix_B[i][j];
        }
    }

   
}



void matrices_add_n(Matrix_t matrix_A,Matrix_t matrix_B){
    assert(matrix_A.width==matrix_B.width && "The dimension of the matrices do not match");
    assert(matrix_A.height==matrix_B.height && "The dimension of the matrices do not match");
    for(int i = 0;i< matrix_A.height;i++){
        for(int j = 0;j< matrix_B.width;j++){
          matrix_A.matrix[i][j] = matrix_A.matrix[i][j]+matrix_B.matrix[i][j];
        }
    }

   
}

void transpose_n(Matrix_t matrix_A,Matrix_t transposed){
  for(int i = 0;i< matrix_A.height;i++){
        for(int j = 0;j< matrix_A.width;j++){
          transposed.matrix[j][i] = matrix_A.matrix[i][j];
          //printf("%f \n",transposed.matrix[j][i] );
        }
    }
  }



void matrices_mul_n(
    Matrix_t matrix_A,
    Matrix_t matrix_B){
    assert(matrix_A.width==matrix_B.width && "The dimension of the matrices do not match");
    assert(matrix_A.height==matrix_B.height && "The dimension of the matrices do not match");
    for(int i = 0;i< matrix_A.height;i++){
        for(int j = 0;j< matrix_A.width;j++){
          matrix_A.matrix[i][j] = matrix_A.matrix[i][j]*matrix_B.matrix[i][j];
        }
    }


}
void matrices_sub_n(
    Matrix_t matrix_A,
    Matrix_t matrix_B){
    assert(matrix_A.width==matrix_B.width && "The dimension of the matrices do not match");
    assert(matrix_A.height==matrix_B.height && "The dimension of the matrices do not match");
    for(int i = 0;i< matrix_A.height;i++){
        for(int j = 0;j< matrix_A.width;j++){
          matrix_A.matrix[i][j] = matrix_A.matrix[i][j]-matrix_B.matrix[i][j];
        }
    }


}
long double matrices_add_all_elements_n(Matrix_t matrix_A){
    long double res = 0;
    for(int i = 0;i< matrix_A.height;i++){
        for(int j = 0;j< matrix_A.width;j++){
         res+=matrix_A.matrix[i][j];
        }
    } 
    return res;
}



void matrices_flatten(
  Image_t matrix, 
  Matrix_t res){
  int counter = 0;
  for(int k = 0;k< matrix.channels;k++){
    for(int i = 0;i< matrix.height;i++){
          for(int j = 0;j< matrix.width;j++){
          res.matrix[counter][0]  = matrix.image[k].matrix[i][j];
          counter++;
          }
      }
  }
 

}



void matrices_reshape_from_flatten_n( 
  Matrix_t matrix, 
  Image_t res){
    int f = 0;
    int g = 0;
    for(int k = 0;k< res.channels;k++){
      for(int i = 0;i< res.height;i++){
        for(int j = 0;j< res.width;j++){
          res.image[k].matrix[i][j]=matrix.matrix[f][g];
          g++;
          if(g>=matrix.width){
             g=0;
             f++; 
          }
        }
      }
    }
  }
void initializeFilter(Filter_t* filt, uint32_t depth_out, uint32_t depth_in, uint32_t f){
  filt->f=f;
  filt->depth_in=depth_in;
  filt->depth_out=depth_out;
  filt->filter = (Matrix_t**)malloc(depth_out*sizeof(Matrix_t*));
  
  for(int i = 0;i<depth_out;i++){
    filt->filter[i] = (Matrix_t*)malloc(depth_out*sizeof(Matrix_t));
    for(int j = 0;j<depth_in;j++){  
      Matrix_t temp;
      matrix_initialize_n(&temp,f,f);
      filt->filter[i][j]=temp;
    } 
  }
}

void image_initialize(
  Image_t* image,
  unsigned int depth,
  unsigned int height,
  unsigned int width){ 
  image->channels=depth;
  image->height=height;
  image->width=width;
  for(int i = 0;i< depth;i++){
      Matrix_t temp;
      
      matrix_initialize_n(&temp,height,width);
      
      image->image[i]=temp;
      
    }

  }
  
void image_initialize_free(
  Image_t image){
  for(int i = 0;i< image.channels;i++){
     matrix_initialize_n_free(image.image[i]);
    }
  }
void initializeFilterFree(Filter_t filt){
  
  
  for(int i = 0;i<filt.depth_out;i++){
    for(int j = 0;j<filt.depth_in;j++){  
      matrix_initialize_n_free(filt.filter[i][j]);
    } 
  }
}

 
void image_to_matrix(//RGB
  bgra_t* src_matrix,
  Image_t* res,
  uint32_t depth, uint32_t height,uint32_t width){
  int c =0;
  res->channels=depth;
  res->height=height;
  res->width=width;
  long double* pixels = (long double *)malloc(res->channels*res->width*res->height*sizeof(long double));
  for(int i  = 0; i<150;i++){
      for(int j  = 0; j<150;j++){
        pixels[c]=src_matrix[i*150+ j].r;
        pixels[c+1]=src_matrix[i*150+ j].g;
        pixels[c+2]=src_matrix[i*150+ j].b;
        c+=3;
        
      }
      
  }


  c=0;
      
  res->image = (Matrix_t*)malloc(res->channels*sizeof(Matrix_t));
  for(int i = 0;i<res->channels;i++){
    Matrix_t temp;
    matrix_initialize_n(&temp,res->height,res->width);
    res->image[i]=temp;
    }
    for (int j = 0; j < res->height; j++){
      for (int k = 0; k < res->width; k++){
        
        res->image[0].matrix[j][k]=pixels[c]/255;
        res->image[1].matrix[j][k]=pixels[c+1]/255;
        res->image[2].matrix[j][k]=pixels[c+2]/255;
        
        c+=3;
      }
    }
  
  free(pixels);

}

void slice_n(
  Image_t matrix_A,
  int start_y,
  int start_x,
  int step,
  Image_t res){

   assert(start_x+step<=matrix_A.width && "Limit out of bounds");
   assert(start_y+step<=matrix_A.height && "Limit out of bounds");

   for(int t =0;t<matrix_A.channels;t++){ 
     int k = 0;
     for(int i = start_y;i<start_y+step;i++){
        int u=0;
        for(int j = start_x;j<start_x+step;j++){
            
            res.image[t].matrix[k][u] = matrix_A.image[t].matrix[i][j]; 
            
            u++;
            
        }    
        k++;
             
     }
   }
   

}






void slice_free_n(
  int start_y,
  int step,
  Image_t res){
  for(int t =0;t<res.channels;t++){ 
    int k=0;
    for(int p = start_y;p<start_y+step;p++){
      free(res.image[t].matrix[k]); 
      k++;
    }
    free(res.image[t].matrix);   
  }

}

void convolution(Image_t image,Filter_t filt,Matrix_t bias,int s, int f, Image_t res,int debugger){
   //assert(image.channels==filt.depth_in && "Dimensions of filter must match dimensions of input image");
  
 
   
   for(int i = 0;i< filt.depth_out;i++){
    int curr_y = 0;
    int out_y = 0;
    
    while(curr_y + f <= image.height){
      int curr_x=0;
      int out_x = 0;
      while(curr_x + f <= image.height){
          long double sum = 0;
         
          Image_t slice_image;   
          slice_image.image = (Matrix_t*)malloc(image.channels*sizeof(Matrix_t));
          
          image_initialize(&slice_image,image.channels,f,f);
          slice_n(image,curr_y,curr_x,f,slice_image);
          for(int j =0;j<image.channels;j++){
            matrices_mul_n(slice_image.image[j],filt.filter[i][j]);
            sum+= matrices_add_all_elements_n(slice_image.image[j]);
                      }
          image_initialize_free(slice_image);
          
          free(slice_image.image);
        
          sum+=bias.matrix[i][0];
        
       
        res.image[i].matrix[out_y][out_x] =sum;
        curr_x+=s;
        out_x+=1;
      }
      curr_y+=s;
      out_y+=1;
    }
   }



}




void maxpool(Image_t image,unsigned int f, unsigned int s, Image_t res){
    ///FILE * fp;
    //fp = fopen ("d.txt", "w+");
    for(unsigned int i = 0;i<image.channels;i++){
        int curr_y=0;
        int out_y=0;
        while(curr_y+f<=image.height){
            int curr_x =0;
            int out_x  =0;
            while(curr_x+f<=image.width){
                Image_t slice_image;   
                slice_image.image = (Matrix_t*)malloc(image.channels*sizeof(Matrix_t));
                
                image_initialize(&slice_image,image.channels,f,f);
                slice_n(image,curr_y,curr_x,f,slice_image);
                
                long double max=slice_image.image[i].matrix[0][0];
                for(int j=0;j<f;j++){
                    for(int k=0;k<f;k++){
                        if(max<slice_image.image[i].matrix[j][k])max=slice_image.image[i].matrix[j][k];
                    }  
                }
                res.image[i].matrix[out_y][out_x]=max;
                //fprintf(fp, "%d\n",(int)max);
                curr_x+=s;
                out_x+=1;

             
                image_initialize_free(slice_image);
                free(slice_image.image);
        

            }
            curr_y+=s;
            out_y+=1;

      
        }

    }
    //
  
  
  //fclose(fp);  
}


void backwardConvolution(Image_t conv_in,Image_t dconv_prev,Filter_t dfilt,Filter_t filt,Matrix_t dbias,int s, int f,Image_t dout){
      //assert(conv_in.channels==filt.depth_in && "Dimensions of filter must match dimensions of input image");

  ////FILE * fp;
    ////fp = fopen ("backConv_.txt", "w+");
   for(int i = 0;i< filt.depth_out;i++){
    int curr_y = 0;
    int out_y = 0;
    
    while(curr_y + f <= conv_in.width){
      int curr_x=0;
      int out_x = 0;
      while(curr_x + f <= conv_in.width){
          
        
        for(int j = 0;j<filt.depth_in;j++){
            for(int k = 0;k<f;k++){
              for(int m = 0;m<f;m++){

                dfilt.filter[i][j].matrix[k][m]+= dconv_prev.image[i].matrix[out_y][out_x]*conv_in.image[j].matrix[curr_y+k][curr_x+m];
              }
            }
          }
        for(int j = 0;j<conv_in.channels;j++){
          for(int k = 0;k<f;k++){
            for(int m = 0;m<f;m++){
              
              //printf("%d %d %d",j,k,m,image_channels);
              dout.image[j].matrix[curr_y+k][curr_x+m]+=(dconv_prev.image[i].matrix[out_y][out_x])*(filt.filter[i][j].matrix[k][m]);
            }
          }
        }
        
       
        curr_x+=s;
        out_x+=1;
      }
      curr_y+=s;
      out_y+=1;
    }
    dbias.matrix[i][0]+=matrices_add_all_elements_n(dconv_prev.image[i]);

    
   }

   //fclose(fp);
  

}



void backwardMaxpool(Image_t dpool,Image_t image,unsigned int f, unsigned int s, Image_t res){
    //FILE * fp;
    //fp = fopen ("backwardMaxpool.txt", "w+");
    for(unsigned int i = 0;i<image.channels;i++){
        int curr_y=0;
        int out_y=0;
        while(curr_y+f<=image.height){
            int curr_x =0;
            int out_x  =0;
            while(curr_x+f<=image.width){
                Image_t slice_image;   
                slice_image.image = (Matrix_t*)malloc(image.channels*sizeof(Matrix_t));
                
                image_initialize(&slice_image,image.channels,f,f);
                slice_n(image,curr_y,curr_x,f,slice_image);
                
                long double max=slice_image.image[i].matrix[0][0];
                
                int a=0;
                int b=0;
                for(int j=0;j<f;j++){
                    for(int k=0;k<f;k++){
                        if(max<slice_image.image[i].matrix[j][k]){
                          max=slice_image.image[i].matrix[j][k];
                          a=j;
                          b=k;
                        }
                    }  
                }
                res.image[i].matrix[curr_y+a][curr_x+b]=dpool.image[i].matrix[out_y][out_x];
                              
               // fprintf(fp, "%d\n",(int)dpool[i][out_y][out_x]);
                
             
                image_initialize_free(slice_image);
                free(slice_image.image);
                curr_x+=s;
                out_x+=1;

            }
            curr_y+=s;
            out_y+=1;

      
        }


    }
  //fclose(fp);  
}

void softmax(Matrix_t x,int size){

  long double sum = 0;
  for(int i=0;i<size;i++)sum+=exp(x.matrix[i][0]);
  for(int i=0;i<size;i++){
    x.matrix[i][0]=exp(x.matrix[i][0])/sum;
    
  }
  
}

long double categoricalCrossEntropy(Matrix_t probs,int** label,int size){ 
  long double res=0;
  for(int i=0;i<size;i++){
    //printf("%f %f\n",probs.matrix[i][0],log(probs.matrix[i][0]));
    res+=(label[i][0]*log(probs.matrix[i][0]))+((1-label[i][0])*log(1-probs.matrix[i][0]));

  }
  return -res/size;


}


int rand_comparison(const void *a, const void *b)
{
    (void)a; (void)b;

    return rand() % 2 ? +1 : -1;
}
void random_batch(int n, int* vektor,int tot_img){
    for(int i=0;i<tot_img;i++)vektor[i]=i;
   
    qsort(vektor, n, sizeof(int), rand_comparison);
    //for(int i=0;i<tot_img;i++)printf("%d\n",vektor[i]);
     
}


void matrices_relu(Image_t image){
  for(int k = 0;k< image.channels;k++){
      for(int i = 0;i< image.height;i++){
          for(int j = 0;j< image.width;j++){
            if (image.image[k].matrix[i][j]<=0)image.image[k].matrix[i][j] =0;
        }
      }
    }
  }
void matrices_drelu(Image_t dimage,Image_t image){
  for(int k = 0;k< dimage.channels;k++){
      for(int i = 0;i< dimage.height;i++){
          for(int j = 0;j< dimage.width;j++){
            if (image.image[k].matrix[i][j]<=0)dimage.image[k].matrix[i][j] =0;
        }
      }
    }
  }
long double apply_conv(Image_t image, int** labels,  int classes,int f,int conv_s, Filter_t filt1 ,Matrix_t bias1, int pool_f,int pool_s,Filter_t filt2,Matrix_t bias2,Matrix_t w3,Matrix_t w4,long double*** grad_w,long double**** grad_f,long double**** grad_b,Matrix_t b3,Matrix_t b4,Filter_t dfilt1,Filter_t dfilt2,Matrix_t dbias1,Matrix_t dbias2,Matrix_t dbias3,Matrix_t dbias4, Matrix_t dw3, Matrix_t dw4, int* predicted,Filter_t filt5,Matrix_t bias5,Filter_t dfilt5,Matrix_t dbias5,Filter_t filt7,Matrix_t bias7,Filter_t dfilt7,Matrix_t dbias7){
    //print_image("image.txt", image);
    Image_t conv1; 
    conv1.image = (Matrix_t*)malloc(filt1.depth_out*sizeof(Matrix_t));
    int conv1_image = (int)((image.height - f)/conv_s)+1;
    image_initialize(&conv1,filt1.depth_out,conv1_image,conv1_image);
    convolution(image,filt1,bias1,conv_s,f,conv1,0);
    matrices_relu(conv1);
    
    
    remove("imagen_gato_c1.txt");
    print_image("imagen_gato_c1.txt",conv1);
   

    //print_image("con1_101.txt", conv1);
    
    unsigned int height_1=(unsigned int)((conv1_image - pool_f)/pool_s)+1;
    unsigned int width_1=(unsigned int)((conv1_image - pool_f)/pool_s)+1;
    Image_t pooled_1; 
    pooled_1.image = (Matrix_t*)malloc(filt1.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled_1,filt1.depth_out,height_1,width_1);
    
    maxpool(conv1, pool_f,  pool_s,pooled_1);

    //print_image("pooled_1_101.txt", pooled_1);
    Image_t conv3; 
    conv3.image = (Matrix_t*)malloc(filt5.depth_out*sizeof(Matrix_t));
    int conv3_image = (int)((pooled_1.height - f)/conv_s)+1;
    image_initialize(&conv3,filt5.depth_out,conv3_image,conv3_image);
    convolution(pooled_1,filt5,bias5, conv_s, f,conv3,0);
    matrices_relu(conv3);
    
    //print_image("conv3_101.txt", conv3);
    unsigned int height_3=(unsigned int)((conv3_image - pool_f)/pool_s)+1;
    unsigned int width_3=(unsigned int)((conv3_image - pool_f)/pool_s)+1;
    Image_t pooled_3; 
    pooled_3.image = (Matrix_t*)malloc(filt5.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled_3,filt5.depth_out,height_3,width_3);
    maxpool(conv3, pool_f,  pool_s,pooled_3);
    
    //print_image("pooled_3_101.txt", pooled_3);
    Image_t conv7; 
    conv7.image = (Matrix_t*)malloc(filt7.depth_out*sizeof(Matrix_t));
    int conv7_image = (int)((pooled_3.height - f)/conv_s)+1;
    image_initialize(&conv7,filt7.depth_out,conv7_image,conv7_image);
    convolution(pooled_3,filt7,bias7, conv_s, f,conv7,0);
    matrices_relu(conv7);
    //print_image("conv7_101.txt", conv7);
    
    unsigned int height_7=(unsigned int)((conv7_image - pool_f)/pool_s)+1;
    unsigned int width_7=(unsigned int)((conv7_image - pool_f)/pool_s)+1;
    Image_t pooled_7; 
    pooled_7.image = (Matrix_t*)malloc(filt7.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled_7,filt7.depth_out,height_7,width_7);
    maxpool(conv7, pool_f,  pool_s,pooled_7);


    //print_image("pooled_7_101.txt", pooled_7);


    Image_t conv2; 
    conv2.image = (Matrix_t*)malloc(filt2.depth_out*sizeof(Matrix_t));
    int conv2_image = (int)((pooled_7.height - f)/conv_s)+1;
    image_initialize(&conv2,filt2.depth_out,conv2_image,conv2_image);
    convolution(pooled_7,filt2,bias2, conv_s, f,conv2,0);
    matrices_relu(conv2);
    //print_image("conv2_101.txt", conv2);

    unsigned int height=(unsigned int)((conv2_image - pool_f)/pool_s)+1;
    unsigned int width=(unsigned int)((conv2_image - pool_f)/pool_s)+1;
    Image_t pooled; 
    pooled.image = (Matrix_t*)malloc(filt2.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled,filt2.depth_out,height,width);
    maxpool(conv2, pool_f,  pool_s,pooled);
   
    //print_image("pooled_101.txt", pooled);
    Matrix_t fc;
    matrix_initialize_n(&fc,filt2.depth_out*width*height,1);
    matrices_flatten( pooled,fc);
    
    Matrix_t z;
    matrix_initialize_n(&z,512,1);
    dot_n(w3,fc,z);
    matrices_add_n(z,b3);  
    //RELU
    for(int i=0;i<512;i++){
      for(int j=0;j<1;j++){
       if(z.matrix[i][j]<=0)z.matrix[i][j]=0;
          
      }    
    }
    
    Matrix_t last;
    
    matrix_initialize_n(&last,2,1);
    dot_n(w4,z,last);

    matrices_add_n(last,b4); 
    softmax(last,2);    
    //print_matrix("last",last);
    //printf("%f %f \n ",last.matrix[0][0],last.matrix[0][1]);
    if (last.matrix[0][0]>last.matrix[1][0]){
      *predicted=0;
    }  else{
      *predicted=1;
    }
    long double loss = categoricalCrossEntropy(last,labels,2);
  /*  ################################################
      ############# Backward Operation ###############
      ################################################ */ 
   
    
    Matrix_t dout;
    
    matrix_initialize_n(&dout,classes,1);
    
  
    for(int i=0;i<classes;i++){
      dout.matrix[i][0]=last.matrix[i][0]-(long double)labels[i][0];
    }
    //printf("%Lf %Lf\n", dout.matrix[0][0],dout.matrix[1][0]);
    Matrix_t z_transpose;
    
    matrix_initialize_n(&z_transpose,z.width,z.height);
    transpose_n(z,z_transpose);
    dot_n(dout,z_transpose,dw4);
    //TODO faltaria hacer db4 pero parecer ser lo mismo que dout
    
    dbias4.matrix[0][0] =dout.matrix[0][0];
    dbias4.matrix[1][0] =dout.matrix[1][0];
    Matrix_t w4_transpose;
    matrix_initialize_n(&w4_transpose,w4.width,w4.height);
        
    
    transpose_n(w4,w4_transpose);
    
    Matrix_t dz;
    matrix_initialize_n(&dz,z.height,z.width);
         
   

    dot_n(w4_transpose,dout,dz);
    
    //print_matrix("z_n",z);
    for(int i=0;i<dz.height;i++){
      for(int j=0;j<dz.width;j++){
       if((z.matrix[i][j]*100000000)<=0)dz.matrix[i][j]=0;
      }    
    }
    
    
    Matrix_t fc_transpose;
    matrix_initialize_n(&fc_transpose,fc.width,fc.height);
    

    transpose_n(fc,fc_transpose);


    

    dot_n(dz,fc_transpose,dw3);
    
    copy_matrix(dbias3,dz);
    
    Matrix_t w3_transpose;
    matrix_initialize_n(&w3_transpose,w3.width,w3.height);
    
    transpose_n(w3,w3_transpose);
    
    Matrix_t dfc;

    matrix_initialize_n(&dfc,fc.height,fc.width);
    
    
    dot_n(w3_transpose,dz,dfc);
  
    Image_t dpool; 
    dpool.image = (Matrix_t*)malloc(filt2.depth_out*sizeof(Matrix_t));
    
    
    image_initialize(&dpool,filt2.depth_out,height,width);
    
    matrices_reshape_from_flatten_n(dfc,dpool);
    //print_image("dpool_101.txt", dpool);
    Image_t dconv2; 
    dconv2.image = (Matrix_t*)malloc(filt2.depth_out*sizeof(Matrix_t));
    image_initialize(&dconv2,filt2.depth_out,conv2_image,conv2_image);
    backwardMaxpool(dpool,conv2,pool_f,pool_s,dconv2);
    
    
    matrices_drelu( dconv2,conv2);
    ////print_image("dconv2_101.txt", dconv2);
    Image_t dpool7; 
    dpool7.image = (Matrix_t*)malloc(filt7.depth_out*sizeof(Matrix_t));
    image_initialize(&dpool7,filt7.depth_out,conv7_image,conv7_image);
    
    backwardConvolution(pooled_7,dconv2,dfilt2,filt2,dbias2,conv_s,f, dpool7);
    //print_filter("dpool7_101.txt", dfilt2);
    Image_t dconv7; 
    dconv7.image = (Matrix_t*)malloc(filt7.depth_out*sizeof(Matrix_t));
    image_initialize(&dconv7,filt7.depth_out,conv7_image,conv7_image);
    backwardMaxpool(dpool7,conv7,pool_f,pool_s,dconv7);
    
    matrices_drelu( dconv7,conv7);
    //print_image("dconv7_101.txt", dconv7);
    Image_t dpool3; 
    dpool3.image = (Matrix_t*)malloc(filt5.depth_out*sizeof(Matrix_t));
    image_initialize(&dpool3,filt5.depth_out,conv3_image,conv3_image);
    backwardConvolution(pooled_3,dconv7,dfilt7,filt7,dbias7,conv_s,f, dpool3);
    //    print_filter("dpool3_101.txt", dfilt7);
    Image_t dconv3; 
    dconv3.image = (Matrix_t*)malloc(filt5.depth_out*sizeof(Matrix_t));
    image_initialize(&dconv3,filt5.depth_out,conv3_image,conv3_image);
    backwardMaxpool(dpool3,conv3,pool_f,pool_s,dconv3);
    
    matrices_drelu( dconv3,conv3);
    //print_image("dconv3_101.txt", dconv3);
    Image_t dpool1; 
    dpool1.image = (Matrix_t*)malloc(filt1.depth_out*sizeof(Matrix_t));
    image_initialize(&dpool1,filt1.depth_out,conv1_image,conv1_image);
    backwardConvolution(pooled_1,dconv3,dfilt5,filt5,dbias5,conv_s,f, dpool1);
    
    //print_filter("dpool1_101.txt", dfilt5);

    Image_t dconv1; 
    dconv1.image = (Matrix_t*)malloc(filt1.depth_out*sizeof(Matrix_t));
    image_initialize(&dconv1,filt1.depth_out,conv1_image,conv1_image);
    backwardMaxpool(dpool1,conv1,pool_f,pool_s,dconv1);
    

    matrices_drelu( dconv1,conv1);
    //print_image("dconv1_101.txt", dconv1);
    Image_t dimage; 
    dimage.image = (Matrix_t*)malloc(filt1.depth_out*sizeof(Matrix_t));
    image_initialize(&dimage,image.channels,image.width,image.width);
    backwardConvolution(image,dconv1,dfilt1,filt1,dbias1,conv_s,f, dimage);
    //print_filter("dimage_101.txt", dfilt1);
    image_initialize_free(conv1);
    image_initialize_free(conv2);
    image_initialize_free(conv3);
    image_initialize_free(conv7);
    image_initialize_free(pooled); 
    image_initialize_free(dpool);
    image_initialize_free(pooled_1);
    image_initialize_free(pooled_3);
    image_initialize_free(pooled_7); 
    image_initialize_free(dpool1);   
    image_initialize_free(dpool3);   
    image_initialize_free(dpool7);  
    image_initialize_free(dconv7); 
    image_initialize_free(dconv3);
    image_initialize_free(dconv2);  
    image_initialize_free(dconv1);  
    image_initialize_free(dimage);  
    matrix_initialize_n_free(fc);
    matrix_initialize_n_free(z);
    matrix_initialize_n_free(last);
    matrix_initialize_n_free(dout);
    matrix_initialize_n_free(z_transpose);
    matrix_initialize_n_free(w4_transpose);
    matrix_initialize_n_free(dz);
    matrix_initialize_n_free(fc_transpose);
    matrix_initialize_n_free(w3_transpose);
    matrix_initialize_n_free(dfc);
      
    return loss;
}
                
long double apply_conv_test(Image_t image, int** labels,  int classes,int f,int conv_s, Filter_t filt1 ,Matrix_t bias1, int pool_f,int pool_s,Filter_t filt2,Matrix_t bias2,Matrix_t w3,Matrix_t w4,Matrix_t b3,Matrix_t b4,int* predicted,Filter_t filt5,Matrix_t bias5,Filter_t filt7,Matrix_t bias7){
    
    //print_image("image.txt", image);
    Image_t conv1; 
    conv1.image = (Matrix_t*)malloc(filt1.depth_out*sizeof(Matrix_t));
    int conv1_image = (int)((image.height - f)/conv_s)+1;
    image_initialize(&conv1,filt1.depth_out,conv1_image,conv1_image);
    convolution(image,filt1,bias1,conv_s,f,conv1,0);
    matrices_relu(conv1);
    //print_image("con1_101.txt", conv1);
    
    unsigned int height_1=(unsigned int)((conv1_image - pool_f)/pool_s)+1;
    unsigned int width_1=(unsigned int)((conv1_image - pool_f)/pool_s)+1;
    Image_t pooled_1; 
    pooled_1.image = (Matrix_t*)malloc(filt1.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled_1,filt1.depth_out,height_1,width_1);
    
    maxpool(conv1, pool_f,  pool_s,pooled_1);

    //print_image("pooled_1_101.txt", pooled_1);
    Image_t conv3; 
    conv3.image = (Matrix_t*)malloc(filt5.depth_out*sizeof(Matrix_t));
    int conv3_image = (int)((pooled_1.height - f)/conv_s)+1;
    image_initialize(&conv3,filt5.depth_out,conv3_image,conv3_image);
    convolution(pooled_1,filt5,bias5, conv_s, f,conv3,0);
    matrices_relu(conv3);
    
    //print_image("conv3_101.txt", conv3);
    unsigned int height_3=(unsigned int)((conv3_image - pool_f)/pool_s)+1;
    unsigned int width_3=(unsigned int)((conv3_image - pool_f)/pool_s)+1;
    Image_t pooled_3; 
    pooled_3.image = (Matrix_t*)malloc(filt5.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled_3,filt5.depth_out,height_3,width_3);
    maxpool(conv3, pool_f,  pool_s,pooled_3);
    
    //print_image("pooled_3_101.txt", pooled_3);
    Image_t conv7; 
    conv7.image = (Matrix_t*)malloc(filt7.depth_out*sizeof(Matrix_t));
    int conv7_image = (int)((pooled_3.height - f)/conv_s)+1;
    image_initialize(&conv7,filt7.depth_out,conv7_image,conv7_image);
    convolution(pooled_3,filt7,bias7, conv_s, f,conv7,0);
    matrices_relu(conv7);
    //print_image("conv7_101.txt", conv7);
    
    unsigned int height_7=(unsigned int)((conv7_image - pool_f)/pool_s)+1;
    unsigned int width_7=(unsigned int)((conv7_image - pool_f)/pool_s)+1;
    Image_t pooled_7; 
    pooled_7.image = (Matrix_t*)malloc(filt7.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled_7,filt7.depth_out,height_7,width_7);
    maxpool(conv7, pool_f,  pool_s,pooled_7);


    //print_image("pooled_7_101.txt", pooled_7);


    Image_t conv2; 
    conv2.image = (Matrix_t*)malloc(filt2.depth_out*sizeof(Matrix_t));
    int conv2_image = (int)((pooled_7.height - f)/conv_s)+1;
    image_initialize(&conv2,filt2.depth_out,conv2_image,conv2_image);
    convolution(pooled_7,filt2,bias2, conv_s, f,conv2,0);
    matrices_relu(conv2);
    //print_image("conv2_101.txt", conv2);

    unsigned int height=(unsigned int)((conv2_image - pool_f)/pool_s)+1;
    unsigned int width=(unsigned int)((conv2_image - pool_f)/pool_s)+1;
    Image_t pooled; 
    pooled.image = (Matrix_t*)malloc(filt2.depth_out*sizeof(Matrix_t));
    image_initialize(&pooled,filt2.depth_out,height,width);
    maxpool(conv2, pool_f,  pool_s,pooled);
   
    //print_image("pooled_101.txt", pooled);
    Matrix_t fc;
    matrix_initialize_n(&fc,filt2.depth_out*width*height,1);
    matrices_flatten( pooled,fc);
    
    Matrix_t z;
    matrix_initialize_n(&z,512,1);
    dot_n(w3,fc,z);
    matrices_add_n(z,b3);  
    //RELU
    for(int i=0;i<512;i++){
      for(int j=0;j<1;j++){
       if(z.matrix[i][j]<=0)z.matrix[i][j]=0;
          
      }    
    }
    
    Matrix_t last;
    
    matrix_initialize_n(&last,2,1);
    dot_n(w4,z,last);

    matrices_add_n(last,b4); 
    softmax(last,2);    
    //print_matrix("last",last);
    //printf("%f %f \n ",last.matrix[0][0],last.matrix[0][1]);
    if (last.matrix[0][0]>last.matrix[1][0]){
      *predicted=0;
    }  else{
      *predicted=1;
    }
    long double loss = categoricalCrossEntropy(last,labels,2);
    image_initialize_free(conv1);
    image_initialize_free(conv2);
    image_initialize_free(conv3);
    image_initialize_free(conv7);
    image_initialize_free(pooled);
    image_initialize_free(pooled_1);
    image_initialize_free(pooled_3);
    image_initialize_free(pooled_7); 
    matrix_initialize_n_free(fc);
    matrix_initialize_n_free(z);
    matrix_initialize_n_free(last);
      
    return loss;
}


void adam_filter(Filter_t f,Filter_t v,Filter_t s, Filter_t d, long double beta1, long double beta2, long double lr,int debugger,int bias_correction_term){
  Filter_t d_copy_2;
  initializeFilter(&d_copy_2,d.depth_out,d.depth_in,d.f);
  copy_filter(d_copy_2,d);
  for(int i =0;i<d_copy_2.depth_out;i++){
    for(int j =0;j<d_copy_2.depth_in;j++){
      matrix_mul_scalar(v.filter[i][j],beta1);//self.beta1*self.m_dw
      matrix_mul_scalar(d_copy_2.filter[i][j],1-beta1);//(1-self.beta1)*dw
      matrices_add_n(v.filter[i][j],d_copy_2.filter[i][j]);//self.beta1*self.m_dw + (1-self.beta1)*dw
      }
  }

  Filter_t d_copy;
  initializeFilter(&d_copy,d.depth_out,d.depth_in,d.f);
  copy_filter(d_copy,d);
  for(int i =0;i<d_copy.depth_out;i++){
    for(int j =0;j<d_copy.depth_in;j++){
      matrix_mul_scalar(s.filter[i][j],beta2);//self.beta2*self.v_dw 
      matrices_mul_n(d_copy.filter[i][j],d_copy.filter[i][j]);//dw**2
      matrix_mul_scalar(d_copy.filter[i][j],1-beta2);//(1-self.beta2)*(dw**2)
      matrices_add_n(s.filter[i][j],d_copy.filter[i][j]);//self.beta2*self.v_dw + (1-self.beta2)*(dw**2)
      }
  }

  Filter_t s_copy;
  initializeFilter(&s_copy,s.depth_out,s.depth_in,s.f);
  Filter_t v_copy;
  initializeFilter(&v_copy,v.depth_out,v.depth_in,v.f);
  copy_filter(s_copy,s);
  copy_filter(v_copy,v);
  for(int i =0;i<s.depth_out;i++){
    for(int j =0;j<s.depth_in;j++){
     
        matrix_div_scalar(s_copy.filter[i][j],1-pow(beta2,bias_correction_term));  //self.v_dw/(1-self.beta2**t)
        matrix_div_scalar(v_copy.filter[i][j],1-pow(beta1,bias_correction_term));  //self.m_dw/(1-self.beta1**t)
      
      }
    }      

  for(int i =0;i<s.depth_out;i++){
    for(int j =0;j<s.depth_in;j++){
      matrix_sqrt(s_copy.filter[i][j],debugger);//np.sqrt(v_dw_corr)
      matrix_add_scalar(s_copy.filter[i][j],0.0000001); //(np.sqrt(v_dw_corr)+self.epsilon) 
      matrix_div(v_copy.filter[i][j],s_copy.filter[i][j]);//(m_dw_corr/(np.sqrt(v_dw_corr)+self.epsilon))
      matrix_mul_scalar(v_copy.filter[i][j],lr);//self.eta*(m_dw_corr/(np.sqrt(v_dw_corr)+self.epsilon))
      matrices_sub_n(f.filter[i][j],v_copy.filter[i][j]);//w - self.eta*(m_dw_corr/(np.sqrt(v_dw_corr)+self.epsilon))
      
      }
    }  
    initializeFilterFree(d_copy_2);
    initializeFilterFree(d_copy);
    initializeFilterFree(s_copy);
    initializeFilterFree(v_copy);
  }

void adam_matrix(Matrix_t w,Matrix_t v,Matrix_t s, Matrix_t d, long double beta1, long double beta2, long double lr,int debugger, int bias_correction_term){
  Matrix_t d_copy;
  matrix_initialize_n(&d_copy,d.height,d.width);
  copy_matrix(d_copy,d);
  Matrix_t d_copy_2;
  matrix_initialize_n(&d_copy_2,d.height,d.width);
  copy_matrix(d_copy_2,d);

  matrix_mul_scalar(v,beta1);//self.beta1*self.m_db
  matrix_mul_scalar(d_copy_2,1-beta1);//(1-self.beta1)*db
  matrices_add_n(v,d_copy_2);//self.beta1*self.m_db + (1-self.beta1)*db  

  matrix_mul_scalar(s,beta2);// self.beta2*self.v_db
  matrices_mul_n(d_copy,d_copy);
  matrix_mul_scalar(d_copy,1-beta2);//(1-self.beta2)*(db)
  matrices_add_n(s,d_copy);//self.beta2*self.v_db + (1-self.beta2)*(db)
  
  Matrix_t s_copy;
  matrix_initialize_n(&s_copy,s.height,s.width);
  copy_matrix(s_copy,s);
  Matrix_t v_copy;
  matrix_initialize_n(&v_copy,v.height,v.width);
  copy_matrix(v_copy,v);
  matrix_div_scalar(s_copy,1-pow(beta2,bias_correction_term));  //self.v_db/(1-self.beta2**t)
  matrix_div_scalar(v_copy,1-pow(beta1,bias_correction_term));  //self.m_db/(1-self.beta1**t)
  
  matrix_sqrt(s_copy,0);
  matrix_add_scalar(s_copy,0.0000001);  
  
  matrix_div(v_copy,s_copy);//(m_db_corr/(np.sqrt(v_db_corr)+self.epsilon))
  matrix_mul_scalar(v_copy,lr);//self.eta*(m_db_corr/(np.sqrt(v_db_corr)+self.epsilon))    
  matrices_sub_n(w,v_copy);
  
  matrix_initialize_n_free(d_copy);
  matrix_initialize_n_free(d_copy_2);
  matrix_initialize_n_free(s_copy);
  matrix_initialize_n_free(v_copy);
}
  


long double optimizer(Image_t* batch_images,int batch, int classes, long double lr,long double beta1, long double beta2, Filter_t* filters,Matrix_t* weights,Matrix_t* biases, Filter_t* filters_grad,Matrix_t* weights_grad,Matrix_t* biases_grad, int itr,int* labels,int conv_s, int pool_s,int f,int pool_f,Filter_t* dfilters,Matrix_t* dbiases,Matrix_t* dweights,int* correct, int bias_correction_term){
  long double cost=0;
  // initialize gradients and momentum,RMS params
  
  
  Filter_t tot_f1;
  
  initializeFilter(&tot_f1,filters[0].depth_out,filters[0].depth_in,filters[0].f);
  
  Matrix_t tot_b1;
  matrix_initialize_n(&tot_b1,biases[0].height,biases[0].width);
  


  Filter_t tot_f2;
  initializeFilter(&tot_f2,filters[1].depth_out,filters[1].depth_in,filters[1].f);
  
  Matrix_t tot_b2;
  matrix_initialize_n(&tot_b2,biases[1].height,biases[1].width);

  Filter_t tot_f5;
  initializeFilter(&tot_f5,filters[2].depth_out,filters[2].depth_in,filters[2].f);
  
  Matrix_t tot_b5;
  matrix_initialize_n(&tot_b5,biases[4].height,biases[4].width);


  Filter_t tot_f7;
  initializeFilter(&tot_f7,filters[3].depth_out,filters[3].depth_in,filters[3].f);
  
  Matrix_t tot_b7;
  matrix_initialize_n(&tot_b7,biases[5].height,biases[5].width);
  

  Matrix_t tot_w3;
  matrix_initialize_n(&tot_w3,weights[0].height,weights[0].width);
  
  Matrix_t tot_w4;
  matrix_initialize_n(&tot_w4,weights[1].height,weights[1].width);

  Matrix_t tot_b3;
  matrix_initialize_n(&tot_b3,biases[2].height,biases[2].width);
  
  Matrix_t tot_b4;
  matrix_initialize_n(&tot_b4,biases[3].height,biases[3].width);
  
  /////
  for(int i =0;i<batch;i++){
    int** y = (int **)malloc( 2*sizeof(int*));
    y[0]=(int *)malloc(sizeof(int));
    y[1]=(int *)malloc(sizeof(int));
    
    if(labels[i]==0){

      y[0][0]=1;

      y[1][0]=0;

    }else{

      y[0][0]=0;

      y[1][0]=1;

    }


    int n_w=2;
    long double*** grad_w = (long double ***)malloc(n_w*sizeof(long double**));
    int n_f=4;
    long double**** grad_f = (long double ****)malloc(n_f*sizeof(long double***));  
      
    int n_b=6;
    long double**** grad_b = (long double ****)malloc(n_b*sizeof(long double***));  
  Filter_t copy_df1;
  initializeFilter(&copy_df1,dfilters[0].depth_out,dfilters[0].depth_in,dfilters[0].f);
  copy_filter(copy_df1,dfilters[0]);
  Matrix_t copy_db1;
  matrix_initialize_n(&copy_db1,dbiases[0].height,dbiases[0].width);
  copy_matrix(copy_db1,dbiases[0]);

  Filter_t copy_df2;
  initializeFilter(&copy_df2,dfilters[1].depth_out,dfilters[1].depth_in,dfilters[1].f);
  copy_filter(copy_df2,dfilters[1]);
  Matrix_t copy_db2;
  matrix_initialize_n(&copy_db2,dbiases[1].height,dbiases[1].width);
  copy_matrix(copy_db2,dbiases[1]);

  Filter_t copy_df5;
  initializeFilter(&copy_df5,dfilters[2].depth_out,dfilters[2].depth_in,dfilters[2].f);
  copy_filter(copy_df5,dfilters[2]);
  Matrix_t copy_db5;
  matrix_initialize_n(&copy_db5,dbiases[4].height,dbiases[4].width);
  copy_matrix(copy_db5,dbiases[4]);

  Filter_t copy_df7;
  initializeFilter(&copy_df7,dfilters[3].depth_out,dfilters[3].depth_in,dfilters[3].f);
  copy_filter(copy_df7,dfilters[3]);
  Matrix_t copy_db7;
  matrix_initialize_n(&copy_db7,dbiases[5].height,dbiases[5].width);
  copy_matrix(copy_db7,dbiases[5]);

  Matrix_t copy_dw3;
  matrix_initialize_n(&copy_dw3,dweights[0].height,dweights[0].width);
  
 
  Matrix_t copy_dw4;
  matrix_initialize_n(&copy_dw4,dweights[1].height,dweights[1].width);

  copy_matrix(copy_dw4,dweights[1]);
  Matrix_t copy_db3;
  matrix_initialize_n(&copy_db3,dbiases[2].height,dbiases[2].width);
  copy_matrix(copy_db3,dbiases[2]);
  Matrix_t copy_db4;
  matrix_initialize_n(&copy_db4,dbiases[3].height,dbiases[3].width);    
  copy_matrix(copy_db4,dbiases[3]);

  int predicted;

  long double loss=apply_conv(batch_images[i], y, classes,f,conv_s,filters[0],biases[0],pool_f,pool_s,filters[1],biases[1],weights[0],weights[1],grad_w,grad_f,grad_b,biases[2],biases[3],copy_df1,copy_df2,copy_db1,copy_db2,copy_db3,copy_db4,copy_dw3,copy_dw4,&predicted,filters[2],biases[4],copy_df5,copy_db5,filters[3],biases[5],copy_df7,copy_db7);
 
  if(predicted==labels[i] )(*correct)++;
    cost+=loss;

    for(int i = 0;i<copy_df1.depth_out;i++){
    for(int j = 0;j<copy_df1.depth_in;j++){  
        
        matrices_add_n(tot_f1.filter[i][j],copy_df1.filter[i][j]);
        
      }
    }  
    for(int i = 0;i<copy_df2.depth_out;i++){
    for(int j = 0;j<copy_df2.depth_in;j++){  
        matrices_add_n(tot_f2.filter[i][j],copy_df2.filter[i][j]);
        
      }
    }
    
    for(int i = 0;i<copy_df5.depth_out;i++){
    for(int j = 0;j<copy_df5.depth_in;j++){  
        matrices_add_n(tot_f5.filter[i][j],copy_df5.filter[i][j]);
        
      }
    }
    for(int i = 0;i<copy_df7.depth_out;i++){
    for(int j = 0;j<copy_df7.depth_in;j++){  
        matrices_add_n(tot_f7.filter[i][j],copy_df7.filter[i][j]);
        
      }
    }
  
    matrices_add_n(tot_w3,copy_dw3);
    matrices_add_n(tot_w4,copy_dw4);
    matrices_add_n(tot_b1,copy_db1);
    matrices_add_n(tot_b2,copy_db2);
    matrices_add_n(tot_b3,copy_db3);
    matrices_add_n(tot_b4,copy_db4);
    matrices_add_n(tot_b5,copy_db5);
    matrices_add_n(tot_b7,copy_db7);
     matrix_initialize_n_free(copy_db1);
  matrix_initialize_n_free(copy_db2);
  matrix_initialize_n_free(copy_db3);
  matrix_initialize_n_free(copy_db4);
  matrix_initialize_n_free(copy_db5);
  matrix_initialize_n_free(copy_db7);
  matrix_initialize_n_free(copy_dw3);
  matrix_initialize_n_free(copy_dw4);
  initializeFilterFree(copy_df1);
  initializeFilterFree(copy_df2);
  initializeFilterFree(copy_df5);
  initializeFilterFree(copy_df7);
  free(y[0]);
  free(y[1]);
  free(y);
  free(grad_b);
  free(grad_w);

  free(grad_f);  
  }


  
  adam_filter(filters[0],filters_grad[0],filters_grad[1], tot_f1,  beta1,  beta2,  lr,0,bias_correction_term);
  adam_matrix(biases[0],biases_grad[0],biases_grad[1], tot_b1,  beta1,  beta2,  lr,0,bias_correction_term);
    
  adam_filter(filters[1],filters_grad[2],filters_grad[3], tot_f2,  beta1,  beta2,  lr,0,bias_correction_term);
  adam_matrix(biases[1],biases_grad[2],biases_grad[3], tot_b2,  beta1,  beta2,  lr,0,bias_correction_term);
  
  adam_filter(filters[2],filters_grad[4],filters_grad[5], tot_f5,  beta1,  beta2,  lr,0,bias_correction_term);
  adam_matrix(biases[4],biases_grad[8],biases_grad[9], tot_b5,  beta1,  beta2,  lr,0,bias_correction_term);
    
  adam_filter(filters[3],filters_grad[6],filters_grad[7], tot_f7,  beta1,  beta2,  lr,0,bias_correction_term);
  adam_matrix(biases[5],biases_grad[10],biases_grad[11], tot_b7,  beta1,  beta2,  lr,0,bias_correction_term);
  
  adam_matrix(weights[0],weights_grad[0],weights_grad[1], tot_w3,  beta1,  beta2,  lr,0,bias_correction_term);
  adam_matrix(biases[2],biases_grad[4],biases_grad[5], tot_b3,  beta1,  beta2,  lr,0,bias_correction_term);
  
  
  adam_matrix(weights[1],weights_grad[2],weights_grad[3], tot_w4,  beta1,  beta2,  lr,0,bias_correction_term);
  adam_matrix(biases[3],biases_grad[6],biases_grad[7], tot_b4,  beta1,  beta2,  lr,0,bias_correction_term);
  
  printf("%Lf\n",cost/batch );
 
  matrix_initialize_n_free(tot_b1);
  matrix_initialize_n_free(tot_b2);
  matrix_initialize_n_free(tot_b3);
  matrix_initialize_n_free(tot_b4);
  matrix_initialize_n_free(tot_b5);
  matrix_initialize_n_free(tot_b7);
  matrix_initialize_n_free(tot_w3);
  matrix_initialize_n_free(tot_w4);

  initializeFilterFree(tot_f1);
  initializeFilterFree(tot_f2);
  initializeFilterFree(tot_f5);
  initializeFilterFree(tot_f7);
  return cost;


}

long double optimizer_test(Image_t* batch_images,int batch, int classes, long double lr,long double beta1, long double beta2, Filter_t* filters,Matrix_t* weights,Matrix_t* biases, int* labels,int conv_s, int pool_s,int f,int pool_f,int* correct){
  long double cost=0;
  // initialize gradients and momentum,RMS params
  
  for(int i =0;i<batch;i++){
    int** y = (int **)malloc( 2*sizeof(int*));
    y[0]=(int *)malloc(sizeof(int));
    y[1]=(int *)malloc(sizeof(int));
    
    if(labels[i]==0){

      y[0][0]=1;

      y[1][0]=0;

    }else{

      y[0][0]=0;

      y[1][0]=1;

    }



  int predicted;
  long double loss= apply_conv_test(batch_images[i], y, classes,f,conv_s,filters[0],biases[0],pool_f,pool_s,filters[1],biases[1],weights[0],weights[1],biases[2],biases[3],&predicted,filters[2],biases[4],filters[3],biases[5]);
  
  if(predicted==labels[i] )(*correct)++;
    cost+=loss;

  
      free(y[0]);
  free(y[1]);
  free(y);
  }

  return cost;


}
void  load_images_from_memory(Image_t* images_train,uint32_t tot_img_train,char* loc ,uint32_t width, uint32_t height, uint32_t image_channels){
  for(int i = 0;i<tot_img_train;i++){
    configuracion_t config;
    config.dst.width = 0;
    config.bits_src = 64;
    config.bits_dst = 64;
    config.es_video = false;
    config.verbose = false;
    config.frames = false;
    config.nombre = false;
    config.cant_iteraciones = 1;
    config.archivo_entrada = NULL;
    config.archivo_entrada_2 = NULL;
    config.carpeta_salida = ".";
    config.extra_archivo_salida = "";
    char buf[33];
    snprintf(buf, 33, loc, i); // puts string into buffer
    //printf("%s\n", buf); // outputs so you can see it
    config.archivo_entrada=buf;
    imagenes_abrir(&config);
    imagenes_flipVertical(&(&config)->src, src_img);
    buffer_info_t info = (&config)->src;
    uint8_t *src =  (uint8_t*)info.bytes;
    bgra_t* src_matrix = (bgra_t*)src;
    Image_t image__; 
    image_to_matrix(src_matrix,&image__,image_channels,height,width);
    images_train[i]= image__;
    free(src_matrix);

    }
}

void fill_filter(Filter_t filt, char* fn){
    FILE *fp = fopen(fn, "r");
    if (fp == NULL)
    {
        printf("Error: could not open file %s", fn);
        
    }

    // reading line by line, max 256 bytes
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH];
   
    for(int i=0;i<filt.depth_out;i++){
      for(int j=0;j<filt.depth_in;j++){
        for(int k=0;k<filt.f;k++){
          for(int l=0;l<filt.f;l++){
            fgets(buffer, MAX_LENGTH, fp);
            filt.filter[i][j].matrix[k][l]=strtod(buffer, NULL);
            //printf("%f %d\n",strtod(buffer, NULL),a++);
          }
        }
      }
    }
  
    // close the file
    fclose(fp);
}

void fill_matrix(Matrix_t matrix, char* fn){
    FILE *fp = fopen(fn, "r");
    if (fp == NULL)
    {
        printf("Error: could not open file %s", fn);
        
    }
    // reading line by line, max 256 bytes
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH];
    for(int i=0;i<matrix.height;i++){
      for(int j=0;j<matrix.width;j++){
        fgets(buffer, MAX_LENGTH, fp);
        matrix.matrix[i][j]=strtod(buffer, NULL);
        //printf("%f %d\n",matrix.matrix[i][j],a++);
      }
    }
  
    // close the file
    fclose(fp);
}
int main( int argc, char** argv ) {

  


  uint32_t test=0;
  //uint32_t number_filter2=64;
  uint32_t conv_s =1;
  uint32_t pool_s =2;
  uint32_t f = 3;
  uint32_t pool_f=2;
  uint32_t image_channels = 3;
  uint32_t image_dim = 150;
  //uint32_t number_filter  = 64 ;
  //uint32_t depth_filter   = 3 ;
  uint32_t classes = 2;
  long double lr = 0.0001;
  long double beta1=0.9;
  long double beta2=0.999;
  if(test==0){  
    uint32_t tot_img_train_cat = 1000;
    Image_t* images_train_cat = (Image_t*)malloc(tot_img_train_cat*sizeof(Image_t));
    
    uint32_t tot_img_train_dog = 1000;
    Image_t* images_train_dog = (Image_t*)malloc(tot_img_train_dog*sizeof(Image_t));
    
    load_images_from_memory( images_train_dog,tot_img_train_dog,"cats_dogs/dogs/%d.bmp",image_dim,image_dim,image_channels);
    load_images_from_memory( images_train_cat,tot_img_train_cat,"cats_dogs/cats/%d.bmp",image_dim,image_dim,image_channels);


    /*
    INICIO PARAMETROS
    */
    Filter_t f1;
    initializeFilter(&f1,32,image_channels,f);
    fill_filter(f1, "f1.txt");
    Matrix_t b1;
    matrix_initialize_n(&b1,32,1);
    //fill_matrix(b1,"weights_1922/b1_trained_aug_v13");
    Filter_t dfilt1;
    initializeFilter(&dfilt1,32,image_channels,f);
    
    Matrix_t dbias1;
    matrix_initialize_n(&dbias1,32,1);
    
    
    Filter_t f2;
    initializeFilter(&f2,128,128,f);
    fill_filter(f2, "f2.txt");
    Matrix_t b2;
    matrix_initialize_n(&b2,128,1);
    //fill_matrix(b2,"weights_1922/b2_trained_aug_v13");
    Filter_t dfilt2;
    initializeFilter(&dfilt2,128,128,f);

    Matrix_t dbias2;
    matrix_initialize_n(&dbias2,128,1);

    Filter_t f5;
    initializeFilter(&f5,64,32,f);
    fill_filter(f5, "f5.txt");
    Matrix_t b5;
    matrix_initialize_n(&b5,64,1);
    //fill_matrix(b5,"weights_1922/b5_trained_aug_v13");
    Filter_t dfilt5;
    initializeFilter(&dfilt5,64,32,f);

    Matrix_t dbias5;
    matrix_initialize_n(&dbias5,64,1);

    Filter_t f7;
    initializeFilter(&f7,128,64,f);
    fill_filter(f7, "f7.txt");

    Matrix_t b7;
    matrix_initialize_n(&b7,128,1);
    //fill_matrix(b7,"weights_1922/b7_trained_aug_v13");
    Filter_t dfilt7;
    initializeFilter(&dfilt7,128,64,f);

    Matrix_t dbias7;
    matrix_initialize_n(&dbias7,128,1);


    
    Matrix_t w3;
    matrix_initialize_n(&w3,512,6272);
    fill_matrix(w3,"w3.txt");
    
    Matrix_t w4;
    matrix_initialize_n(&w4,2,512);
    fill_matrix(w4,"w4.txt");
    Matrix_t dw3;
    matrix_initialize_n(&dw3,512,6272);
    
    Matrix_t dw4;
    matrix_initialize_n(&dw4,2,512);
    
    for(int i = 0;i<w3.height;i++){
      for(int j = 0;j<w3.width;j++){
        dw3.matrix[i][j]=1;
      }
    } 
    
    
    for(int i = 0;i<w4.height;i++){
      for(int j = 0;j<w4.width;j++){
        dw4.matrix[i][j]=1;
      }
    } 


    Matrix_t b3;
    matrix_initialize_n(&b3,512,1);
    //fill_matrix(b3,"weights_1922/b3_trained_aug_v13");

    Matrix_t db3;
    matrix_initialize_n(&db3,512,1);
    

    Matrix_t b4;
    matrix_initialize_n(&b4,2,1);
    //fill_matrix(b4,"weights_1922/b4_trained_aug_v13");


    Matrix_t db4;
    matrix_initialize_n(&db4,2,1);
    
    
     /*Armo vectores de parametros*/  
    
    Matrix_t* biases=(Matrix_t*)malloc( 6*sizeof(Matrix_t));
    Matrix_t* weights=(Matrix_t*)malloc( 2*sizeof(Matrix_t));
    Filter_t* filters=(Filter_t*)malloc( 4*sizeof(Filter_t));
    Matrix_t* dbiases=(Matrix_t*)malloc( 6*sizeof(Matrix_t));
    Matrix_t* dweights=(Matrix_t*)malloc( 2*sizeof(Matrix_t));
    
    Filter_t* dfilters=(Filter_t*)malloc( 4*sizeof(Filter_t));
    
    biases[0]=b1;
    biases[1]=b2;
    biases[2]=b3;
    biases[3]=b4;
    biases[4]=b5;
    biases[5]=b7;
    
    weights[0]=w3;
    weights[1]=w4;

    filters[0]=f1;
    filters[1]=f2;
    filters[2]=f5;
    filters[3]=f7;

    dbiases[0]=dbias1;
    dbiases[1]=dbias2;
    dbiases[2]=db3;
    dbiases[3]=db4;
    dbiases[4]=dbias5;
    dbiases[5]=dbias7;

    dfilters[0]=dfilt1;
    dfilters[1]=dfilt2;
    dfilters[2]=dfilt5;
    dfilters[3]=dfilt7;
    dweights[0]=dw3;
    dweights[1]=dw4;
    
    
    Matrix_t* biases_grad=(Matrix_t*)malloc( 2*6*sizeof(Matrix_t));
    Matrix_t* weights_grad=(Matrix_t*)malloc( 2*2*sizeof(Matrix_t));
    Filter_t* filters_grad=(Filter_t*)malloc( 2*4*sizeof(Filter_t));
    
    Filter_t v1;
    initializeFilter(&v1,32,image_channels,f);

    Filter_t s1;
    initializeFilter(&s1,32,image_channels,f);

    filters_grad[0]=v1;
    filters_grad[1]=s1;
    
   
    
    Filter_t v2;
    initializeFilter(&v2,128,128,f);

    Filter_t s2;
    initializeFilter(&s2,128,128,f);

    filters_grad[2]=v2;
    filters_grad[3]=s2;

    Filter_t v5;
    initializeFilter(&v5,64,32,f);

    Filter_t s5;
    initializeFilter(&s5,64,32,f);

    filters_grad[4]=v5;
    filters_grad[5]=s5;
    
    Filter_t v7;
    initializeFilter(&v7,128,64,f);

    Filter_t s7;
    initializeFilter(&s7,128,64,f);

    filters_grad[6]=v7;
    filters_grad[7]=s7;



    Matrix_t v3;
    matrix_initialize_n(&v3,512,6272);

    Matrix_t s3;
    matrix_initialize_n(&s3,512,6272);
    
    weights_grad[0]=v3;
    weights_grad[1]=s3;


    Matrix_t v4;
    matrix_initialize_n(&v4,2,512);

    Matrix_t s4;
    matrix_initialize_n(&s4,2,512);

    weights_grad[2]=v4;
    weights_grad[3]=s4;
    
    Matrix_t bv1;
    matrix_initialize_n(&bv1,32,1);

    Matrix_t bs1;
    matrix_initialize_n(&bs1,32,1);
    
    biases_grad[0]=bv1;
    biases_grad[1]=bs1;
    
    Matrix_t bv2;
    matrix_initialize_n(&bv2,128,1);

    Matrix_t bs2;
    matrix_initialize_n(&bs2,128,1);
    biases_grad[2]=bv2;
    biases_grad[3]=bs2;
    
    
    Matrix_t bv3;
    matrix_initialize_n(&bv3,512,1);
    
    Matrix_t bs3;
    matrix_initialize_n(&bs3,512,1);
    
    biases_grad[4]=bv3;
    biases_grad[5]=bs3;

    Matrix_t bv4;
    matrix_initialize_n(&bv4,2,1);
    
    Matrix_t bs4;
    matrix_initialize_n(&bs4,2,1);
    
    biases_grad[6]=bv4;
    biases_grad[7]=bs4;

    Matrix_t bv5;
    matrix_initialize_n(&bv5,64,1);

    Matrix_t bs5;
    matrix_initialize_n(&bs5,64,1);
    biases_grad[8]=bv5;
    biases_grad[9]=bs5;

    Matrix_t bv7;
    matrix_initialize_n(&bv7,128,1);

    Matrix_t bs7;
    matrix_initialize_n(&bs7,128,1);
    biases_grad[10]=bv7;
    biases_grad[11]=bs7;



    
      
    /*
    FIN PARAMETROS 
    */
    printf("%s %d\n","Total de imagenes cargadas:",(tot_img_train_cat+tot_img_train_dog));
    int itr;
    unsigned int epochs=100;
    unsigned int batch = 20;
    unsigned int batch_iterations= (tot_img_train_cat+tot_img_train_dog)/batch;
    //long double cost[epochs*batch];
    Image_t* batch_images = (Image_t*)malloc(batch*sizeof(Image_t));
    int file_update=140;
    int bias_correction_term=1;
    for(int i=0;i<epochs;i++){
      printf("%s %d\n","Epcoch: ",i);
      epch_counter++;
      batch_counter=0;
      int vektor[tot_img_train_cat+tot_img_train_dog];
      int labels[batch];
      int contador=0;    
      random_batch(tot_img_train_cat+tot_img_train_dog, vektor,tot_img_train_cat+tot_img_train_dog);
      int correct_total=0;
      for(int j=0;j<batch_iterations;j++){
        printf("%s %d\n","Batch: ",j);
        batch_counter++;
        for(int k=0;k<batch;k++){
          //Cats from 0 to 999. Dogs from 1000 to 1999
          if(vektor[contador]>=tot_img_train_dog){
            batch_images[k]=images_train_dog[vektor[contador]-tot_img_train_dog];
            labels[k]=1;
          }
          if(vektor[contador]<tot_img_train_cat){
            batch_images[k]=images_train_cat[vektor[contador]];
            labels[k]=0;
          }
          contador++;
          itr=j*k;//Aca tengo que guardar el costo
          
        }

        int correct=0;
        long double loss_tot =optimizer(batch_images,batch, classes, lr,beta1, beta2, filters,weights,biases, filters_grad,weights_grad,biases_grad,itr,labels, conv_s, pool_s,f,pool_f,dfilters,dbiases,dweights,&correct,bias_correction_term);
        correct_total+=correct;
        //cost[i*j]=loss_tot;
         bias_correction_term++;
    
      
      }
     
    printf("%s %Lf\n","Acc",(long double)correct_total/(tot_img_train_cat+tot_img_train_dog));
    char f1_file[33];
    char f2_file[33];
    char f5_file[33];
    char f7_file[33];
    char b1_file[33];
    char b2_file[33];
    char b3_file[33];
    char b4_file[33];
    char b5_file[33];
    char b7_file[33];
    char w3_file[33];
    char w4_file[33];
    snprintf(f1_file, 33, "weights_1922/f1_trained_aug_v%d", file_update);
    snprintf(f2_file, 33, "weights_1922/f2_trained_aug_v%d", file_update);
    snprintf(f5_file, 33, "weights_1922/f5_trained_aug_v%d", file_update);
    snprintf(f7_file, 33, "weights_1922/f7_trained_aug_v%d", file_update);
    snprintf(b1_file, 33, "weights_1922/b1_trained_aug_v%d", file_update);
    snprintf(b2_file, 33, "weights_1922/b2_trained_aug_v%d", file_update);
    snprintf(b3_file, 33, "weights_1922/b3_trained_aug_v%d", file_update);
    snprintf(b4_file, 33, "weights_1922/b4_trained_aug_v%d", file_update);
    snprintf(b5_file, 33, "weights_1922/b5_trained_aug_v%d", file_update);
    snprintf(b7_file, 33, "weights_1922/b7_trained_aug_v%d", file_update);
    snprintf(w3_file, 33, "weights_1922/w3_trained_aug_v%d", file_update);
    snprintf(w4_file, 33, "weights_1922/w4_trained_aug_v%d", file_update);
    print_filter(f1_file,f1);
    print_filter(f2_file,f2);
    print_filter(f5_file,f5);
    print_filter(f7_file,f7);
    print_matrix(b1_file,b1);
    print_matrix(b2_file,b2);
    print_matrix(b3_file,b3);
    print_matrix(b4_file,b4);
    print_matrix(b5_file,b5);
    print_matrix(b7_file,b7);
    print_matrix(w3_file,w3);
    print_matrix(w4_file,w4);
    file_update++;

    }


    free(batch_images);  
  }else{
 for(int t=0;t<50;t++){ 
    
      uint16_t tot_img_test_cat = 500;
      Image_t* images_test_cat = (Image_t*)malloc(tot_img_test_cat*sizeof(Image_t));
      
      uint16_t tot_img_test_dog = 500;
      Image_t* images_test_dog = (Image_t*)malloc(tot_img_test_dog*sizeof(Image_t));
      //"test/cats_bmp/%d.bmp"
      //"test/dogs_bmp/%d.bmp"
      load_images_from_memory( images_test_dog,tot_img_test_dog,"test/dogs_bmp/%d.bmp",image_dim,image_dim,image_channels);
      load_images_from_memory( images_test_cat,tot_img_test_cat,"test/cats_bmp/%d.bmp",image_dim,image_dim,image_channels);
     /*
      INICIO PARAMETROS
      */
          /*
    INICIO PARAMETROS
    */

    char f1_file[33];
    snprintf(f1_file, 33, "weights_1922/f1_trained_aug_v%d", t);
    Filter_t f1;
    initializeFilter(&f1,32,image_channels,f);
    fill_filter(f1, f1_file);

    char b1_file[33];
    snprintf(b1_file, 33, "weights_1922/b1_trained_aug_v%d", t);
    Matrix_t b1;
    matrix_initialize_n(&b1,32,1);
    fill_matrix(b1,b1_file);
    
    
    char f2_file[33];
    snprintf(f2_file, 33, "weights_1922/f2_trained_aug_v%d", t);
    Filter_t f2;
    initializeFilter(&f2,128,128,f);
    fill_filter(f2, f2_file);
    
    char b2_file[33];
    snprintf(b2_file, 33, "weights_1922/b2_trained_aug_v%d", t);
    Matrix_t b2;
    matrix_initialize_n(&b2,128,1);
    fill_matrix(b2,b2_file);
    
    char f5_file[33];
    snprintf(f5_file, 33, "weights_1922/f5_trained_aug_v%d", t);
    Filter_t f5;
    initializeFilter(&f5,64,32,f);
    fill_filter(f5, f5_file);
    
    char b5_file[33];
    snprintf(b5_file, 33, "weights_1922/b5_trained_aug_v%d", t);
    
    Matrix_t b5;
    matrix_initialize_n(&b5,64,1);
    fill_matrix(b5,b5_file);
   
    char f7_file[33];
    snprintf(f7_file, 33, "weights_1922/f7_trained_aug_v%d", t);
    
    Filter_t f7;
    initializeFilter(&f7,128,64,f);
    fill_filter(f7, f7_file);

    char b7_file[33];
    snprintf(b7_file, 33, "weights_1922/b7_trained_aug_v%d", t);
    
    Matrix_t b7;
    matrix_initialize_n(&b7,128,1);
    fill_matrix(b7,b7_file);
  

    char w3_file[33];
    snprintf(w3_file, 33, "weights_1922/w3_trained_aug_v%d", t);
    
    Matrix_t w3;
    matrix_initialize_n(&w3,512,6272);
    fill_matrix(w3,w3_file);
    
    char w4_file[33];
    snprintf(w4_file, 33, "weights_1922/w4_trained_aug_v%d", t);
    Matrix_t w4;
    matrix_initialize_n(&w4,2,512);
    fill_matrix(w4,w4_file);
 

    char b3_file[33];
    snprintf(b3_file, 33, "weights_1922/b3_trained_aug_v%d", t);
    Matrix_t b3;
    matrix_initialize_n(&b3,512,1);
    fill_matrix(b3,b3_file);

    
    char b4_file[33];
    snprintf(b4_file, 33, "weights_1922/b4_trained_aug_v%d", t);
    Matrix_t b4;
    matrix_initialize_n(&b4,2,1);
    fill_matrix(b4,b4_file);


    
     /*Armo vectores de parametros*/  
    
    Matrix_t* biases=(Matrix_t*)malloc( 6*sizeof(Matrix_t));
    Matrix_t* weights=(Matrix_t*)malloc( 2*sizeof(Matrix_t));
    Filter_t* filters=(Filter_t*)malloc( 4*sizeof(Filter_t));
   
    biases[0]=b1;
    biases[1]=b2;
    biases[2]=b3;
    biases[3]=b4;
    biases[4]=b5;
    biases[5]=b7;
    
    weights[0]=w3;
    weights[1]=w4;

    filters[0]=f1;
    filters[1]=f2;
    filters[2]=f5;
    filters[3]=f7;



      
        
      /*
      FIN PARAMETROS 
      */
      
      unsigned int batch = (tot_img_test_cat+tot_img_test_dog);
      Image_t* batch_images = (Image_t*)malloc(batch*sizeof(Image_t));
      int vektor[tot_img_test_cat+tot_img_test_dog];
      int labels[batch];
      int contador=0;    
      random_batch(tot_img_test_cat+tot_img_test_dog, vektor,tot_img_test_cat+tot_img_test_dog);
      int correct_total=0;
      for(int k=0;k<batch;k++){
        //Cats from 0 to 999. Dogs from 1000 to 1999
        if(vektor[contador]>=tot_img_test_dog){
          batch_images[k]=images_test_dog[vektor[contador]-tot_img_test_dog];
          labels[k]=1;
        }
        if(vektor[contador]<tot_img_test_cat){
          batch_images[k]=images_test_cat[vektor[contador]];
          labels[k]=0;
        }
        contador++;
        
      }
      int correct=0;
      long double loss_tot =optimizer_test(batch_images,batch, classes, lr,beta1, beta2, filters,weights,biases, labels, conv_s, pool_s,f,pool_f,&correct);
      correct_total+=correct;
      printf("%s %Lf\n","Acc on test",(long double)correct_total/(tot_img_test_cat+tot_img_test_dog));
      


      free(batch_images);  
    }
  }
  return 0;
}
// gcc -o arq1922 arquitectura1922.c -lm -O3
// ./ejec aa -i asm name.bmp


//https://stackoverflow.com/questions/12747731/ways-to-create-dynamic-matrix-in-c







