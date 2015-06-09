// Faire la rotation de Yaroslavsky

#include "iio.h"
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <fftw3.h>


// déclaration des fonctions
int trans_ligne(float *img,int w,int h,float *vect_trans,int l) ;
int trans_colonne(float *img,int w,int h,float *vect_trans,int l) ;
 
void fourier1dForward(float* in,float* reOut,float* imOut,unsigned int largeur) ;
void fourier1dBackward(float* reIn,float* imIn,float* out,unsigned int largeur) ;


// Définition de la valeur absolue
float absolute(float value) {
  if (value < 0) {
    return -value;
  }
  else {
    return value;  
  }
}


// Code de la Rotation de Yaroslavsky :
int yaroslavsky(float *img,float *img_f,int x,int z,int w_f,int h_f,double *a){
//Arguments : yaroslavsky(image initiale,image rotatée,dimension1 de l'image initiale,dimension2 de l'image initiale,dimension1 de l'image finale,dimension2 de l'image finale,

// On va plonger l'image dans un autre image de dimension le double du max des dimension de l'image finale (on prende le max pour bien traiter les images rectangulaires).
int dimax ;
if (x<=z) {dimax=z;}else{dimax=x;}


float *img_aux=malloc(2*dimax*2*dimax*sizeof(float)); 
float *img_aux_dimax=malloc(2*dimax*2*dimax*sizeof(float));
float ligne_trans[2*dimax] ; //vecteur pour le deuxième shear
float one_col_trans[2*dimax] ;//vecteur pour le premier shear
float two_col_trans[2*dimax] ; //vecteur pour le troisième shear
float pi;
double cosangle ; // cosinus de l'angle de la rotation
double sinangle ; // sinus de l'angle de la rotation
float tandemiangle ; // tangente de la moité de l'angle de la rotation
float cosangleaux ;
pi=3.14159265359 ;



for(int p=0;p<3;p++){

	float meanimg; // moyenne de l'image : permet de limiter l'effet gibbs
meanimg=0.0;
        for (int i=0;i<x;i++){
		for (int j=0;j<z;j++){
                meanimg=meanimg+img[(i+j*x)*3+p] ;
}
}   
	

	for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
			img_aux[i+(2*dimax)*(j)]=meanimg/(x*z*1.0) ;
		}
	}


	

     for (int i=0;i<x;i++){
		for (int j=0;j<z;j++){
		 	img_aux[dimax/2+(dimax-x)/2+i+(2*dimax)*(dimax/2+(dimax-z)/2+j)]=img[(i+j*x)*3+p] ;
		}
	}



// conditions aux bords périodisé

int i_sym,j_sym;
	
		for(int i=dimax-x;i<dimax+x;i++){
			for(int j=dimax-z;j<dimax+z;j++){
                i_sym = i-dimax+x/2;
                while(i_sym<0 || i_sym>x-1){i_sym = (i_sym<0) ? i_sym+x : i_sym-x;}
                j_sym = j-dimax+z/2;
                while(j_sym<0 || j_sym>z-1){j_sym = (j_sym<0) ? j_sym+z : j_sym-z;}
				img_aux[(i+2*dimax*j)]=img[(i_sym+j_sym*x)*3+p];
			}
		}
	







cosangle=a[0] ;
sinangle=a[1] ;
// début : se ramener à un angle entre -pi/4 et pi/4
	
// si l'angle est sup à pi/2 en valeur absolue.

if(cosangle<0)

{  for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		 img_aux_dimax[i+2*dimax*j]=img_aux[(2*dimax-i-1)+2*dimax*(2*dimax-j-1)];
                 
		}
		}
  for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		img_aux[i+2*dimax*j]=img_aux_dimax[i+2*dimax*j];
                 
		}
		}
  cosangle=-cosangle ;
  sinangle=-sinangle ;
  tandemiangle=sinangle/(1.0+cosangle) ;
}

 tandemiangle=sinangle/(1.0+cosangle) ;

// si l'angle est supérieur à pi/4 :
if(absolute(sinangle)>absolute(cosangle))
//if(0==0)
{ if (sinangle>=0)
//if(0==0)
{ for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		 img_aux_dimax[i+2*dimax*j]=img_aux[j+2*dimax*(2*dimax-i-1)];
                 
		}
		}
  for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		img_aux[i+2*dimax*j]=img_aux_dimax[i+2*dimax*j];
                 
		}
		}
  cosangleaux=cosangle ;
  cosangle=sinangle ;
  sinangle=-cosangleaux ;
  tandemiangle=sinangle/(1.0+cosangle) ;
}
else if (sinangle<0)

{
for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		 img_aux_dimax[i+2*dimax*j]=img_aux[(2*dimax-j-1)+2*dimax*(i)];
                
		}
		}
  for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		img_aux[i+2*dimax*j]=img_aux_dimax[i+2*dimax*j];
                 
		}
		}
  cosangleaux=cosangle ;
  cosangle=-sinangle ;
  sinangle=cosangleaux ;
  tandemiangle=sinangle/(1.0+cosangle) ;
}
}
 // fin si l'angle sup à pi/4

// début : se ramener à un angle entre -pi/4 et pi/4




// Première translation
	for (int i=0;i<2*dimax;i++){
	   one_col_trans[i]=i*tandemiangle;

	}


	trans_colonne(img_aux,2*dimax,2*dimax,one_col_trans,2*dimax) ;

        for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		 img_aux_dimax[i+2*dimax*j]=img_aux[i+2*dimax*good_modulus(j+one_col_trans[dimax],2*dimax)];
                 
		}
		}
  for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		img_aux[i+2*dimax*j]=img_aux_dimax[i+2*dimax*j];
                 
		}
		}



// seconde translation
	for (int i=0;i<2*dimax;i++){
	   ligne_trans[i]=(-i*sinangle) ;

	}


	trans_ligne(img_aux,2*dimax,2*dimax,ligne_trans,2*dimax) ;
    

	for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		 img_aux_dimax[i+2*dimax*j]=img_aux[good_modulus(i+ligne_trans[dimax],2*dimax)+2*dimax*j];
                 
		}
		}
  for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		img_aux[i+2*dimax*j]=img_aux_dimax[i+2*dimax*j];
                 
		}
		}


// troisième translation
	for (int i=0;i<2*dimax;i++){
	   two_col_trans[i]=i*tandemiangle ;

	}


	trans_colonne(img_aux,2*dimax,2*dimax,two_col_trans,2*dimax) ;

    for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		 img_aux_dimax[i+2*dimax*j]=img_aux[i+2*dimax*good_modulus(j+two_col_trans[dimax],2*dimax)];
                 
		}
		}
  for (int i=0;i<2*dimax;i++){
		for (int j=0;j<2*dimax;j++){
		img_aux[i+2*dimax*j]=img_aux_dimax[i+2*dimax*j];
                 
		}
		}




  

     for (int i=0;i<w_f;i++){
		for (int j=0;j<h_f;j++){
		 img_f[(i+j*w_f)*3+p]=img_aux[(i+dimax-w_f/2)+2*dimax*(j+dimax-h_f/2)];
                 
		}
	}




}
	return 0;
}





// Translation ligne par Fourier :

int trans_ligne(float *img,int w,int h,float *vect_trans,int l){

	        float pi=3.14159265359 ;
		float *refftimg=malloc(w*h*sizeof(float));
		float *imfftimg=malloc(w*h*sizeof(float));
		float *refftrans=malloc(w*h*sizeof(float));
		float ligneimg[w];
                float lignerefftimg[w];
		float ligneimfftimg[w];
		

 		for (int j=0;j<h;j++){
                  for (int i=0;i<w;i++){
		      ligneimg[i]=img[i+w*j] ;
			}
                  fourier1dForward(ligneimg,lignerefftimg,ligneimfftimg,w) ;
                   for (int i=0;i<w;i++){
		      refftimg[i+w*j]=lignerefftimg[i] ;
                      imfftimg[i+w*j]=ligneimfftimg[i] ;
			}

		
		   }
	       


		for (int i=0;i<w;i++){
		   for (int j=0;j<h;j++){
			int ii=i-w/2 ;
                 refftrans[i+j*w]=cos(-2*pi*ii*(vect_trans[j]/((float) w)))*refftimg[i+j*w]-sin(-2*pi*ii*(vect_trans[j]/((float) w)))*imfftimg[i+j*w];
                 imfftimg[i+j*w]=sin(-2*pi*ii*vect_trans[j]/((float) w))*refftimg[i+j*w]+cos(-2*pi*ii*(vect_trans[j]/((float) w)))*imfftimg[i+j*w];		 	
                    
 

		}
		  }

	for (int j=0;j<h;j++){
		for (int i=0;i<w;i++){
		      lignerefftimg[i]=refftrans[i+w*j] ;
                      ligneimfftimg[i]=imfftimg[i+w*j] ;

			}
                  
                  fourier1dBackward(lignerefftimg,ligneimfftimg,ligneimg,w) ;
                   
		for (int i=0;i<w;i++){
		      img[i+w*j]=ligneimg[i] ;
			}
		
		   }
              


	

	return 0;
 }


// Translation colonnes par Fourier:

int trans_colonne(float *img,int w,int h,float *vect_trans,int l){

	        float pi=3.14159265359 ;
	        float *refftimg=malloc(w*h*sizeof(float));
		float *imfftimg=malloc(w*h*sizeof(float));
		float *refftrans=malloc(w*h*sizeof(float));
                float ligneimg[h];
                float lignerefftimg[h];
		float ligneimfftimg[h];

                for (int i=0;i<w;i++){
                  for (int j=0;j<h;j++){
		      ligneimg[j]=img[i+w*j] ;
			}
                  fourier1dForward(ligneimg,lignerefftimg,ligneimfftimg,h) ;
                   for (int j=0;j<h;j++){
		      refftimg[i+w*j]=lignerefftimg[j] ;
                      imfftimg[i+w*j]=ligneimfftimg[j] ;
			}

		
		   }
	       
	
		


		for (int i=0;i<w;i++){
		   for (int j=0;j<h;j++){
		int jj=j-h/2 ;
                 refftrans[i+j*w]=cos(-2.0*pi*jj*vect_trans[i]/((float)h))*refftimg[i+j*w]-sin(-2.0*pi*jj*vect_trans[i]/((float)h))*imfftimg[i+j*w];		 	
                 imfftimg[i+j*w]=sin(-2.0*pi*jj*vect_trans[i]/((float)h))*refftimg[i+j*w]+cos(-2.0*pi*jj*vect_trans[i]/((float)h))*imfftimg[i+j*w];   
                 

		}
		  }

               for (int i=0;i<w;i++){
		for (int j=0;j<h;j++){
		      lignerefftimg[j]=refftrans[i+w*j] ;
                      ligneimfftimg[j]=imfftimg[i+w*j] ;  
			}
                  
                  fourier1dBackward(lignerefftimg,ligneimfftimg,ligneimg,h) ;
                   
		for (int j=0;j<h;j++){
		      img[i+w*j]=ligneimg[j] ;
			}
		
		   }


	

	return 0;
 }



		
 
void fourier1dForward(float* in,
                    float* reOut,
                    float* imOut,
                    unsigned int largeur)
{
   double* spatial_repr;
   fftw_complex* frequency_repr;
   unsigned int i;
   fftw_plan plan;
   int x;

     spatial_repr= malloc(sizeof(double)*(largeur));
     frequency_repr= fftw_malloc(sizeof(fftw_complex)*largeur);
   
 
  

 
   for(i=0;i<largeur;i++)
   {
     
      spatial_repr[i] = (double)(in[i]);
     
   }
 
   /*on calcule le plan d'exécution*/
   plan=fftw_plan_dft_r2c_1d(largeur, spatial_repr, frequency_repr, FFTW_ESTIMATE);
 
   /*on calcule la transformée*/
   fftw_execute(plan);
 

      for(i=0;i<largeur;i++)
      {
	        /*on recentre l'image*/
	      x=good_modulus(i+largeur/2,largeur);
	      
          reOut[x]=frequency_repr[i][0];
          imOut[x]=frequency_repr[i][1];
      }
 
   fftw_destroy_plan(plan);
   fftw_free(spatial_repr);
   fftw_free(frequency_repr);
 
}





void fourier1dBackward(float* reIn,
                     float* imIn,
                     float* out,
                     unsigned int largeur)
{
   double* spatial_repr;
   fftw_complex* frequency_repr;
   unsigned int i;
   int x;
   fftw_plan plan;
 
   //spatial_repr= fftw_malloc(sizeof(fftw_complex)*largeur);
   spatial_repr= malloc(sizeof(double)*largeur);
   frequency_repr= fftw_malloc(sizeof(fftw_complex)*largeur);
 
  
      for(i=0;i<largeur;i++)
      {
          /*on décentre*/
	      x=i;
	  
	      x=good_modulus(i+largeur/2,largeur);

	      frequency_repr[i][0]=reIn[x];
	      frequency_repr[i][1]=imIn[x];
      }
 
  plan=fftw_plan_dft_c2r_1d(largeur, frequency_repr, spatial_repr, FFTW_ESTIMATE);
 
  fftw_execute(plan);
 
   /*on retranscrit l'image complexe en image réelle, sans oublier de diviser par largeur*hauteur*/
   for(i=0;i<largeur;i++)
   {
      
      out[i]=(float)((spatial_repr[i]))/((float)largeur) ;

   }
 
   fftw_destroy_plan(plan);
   fftw_free(spatial_repr);
   fftw_free(frequency_repr);
}
