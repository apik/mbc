#ifndef __CGAMMA_H__
#define __CGAMMA_H__
/*
  
  Conversion of CERNLIB fortran routines 
  WGAMMA and  WPSIPG from Fortran to C
  
  by A.Pikelner
  
*/
#include <complex.h>
#include <stdlib.h>
#include <math.h>

double complex wgamma(double complex Z)
{
  static double C[16];
  static double Z1 = 1.;
  static double HF = 1./2.;
  static double C1 = 2.50662827463100050;
  
  C[ 0] = 41.624436916439068;
  C[ 1] =-51.224241022374774;
  C[ 2] = 11.338755813488977;
  C[ 3] = -0.747732687772388;
  C[ 4] =  0.008782877493061;
  C[ 5] = -0.000001899030264;
  C[ 6] =  0.000000001946335;
  C[ 7] = -0.000000000199345;
  C[ 8] =  0.000000000008433;
  C[ 9] =  0.000000000001486;
  C[10] = -0.000000000000806;
  C[11] =  0.000000000000293;
  C[12] = -0.000000000000102;
  C[13] =  0.000000000000037;
  C[14] = -0.000000000000014;
  C[15] =  0.000000000000006;
  
  double complex U,V,F,H,S;
    
  U=Z;
  double X = creal(U);
  
  if(cimag(U) == 0 && -fabs(X) == (int)X) 
    {
      F=0;
      H=0;
      printf("ARGUMENT EQUALS NON-POSITIVE INTEGER = %f", creal(Z));
      exit(1);
    }
  else
    {
      if(creal(X) >= 1)
        {
          F=1.;
          V=U;
        }
      else if(creal(X) >= 0) 
        {
          F=1./U;
          V=1+U;
        }
      else
        {
          F=1;
          V=1-U;
        }
      H=1.;
      S=C[0];

      int K;
      for( K = 1; K <= 15; K++)
        {
          H=((V-K)/(V+(K-1)))*H;
          S=S+C[K]*H;
        }
        
      H=V+(4+HF);
      H=C1*cexp((V-HF)*clog(H)-H)*S;
      if(creal(X) < 0) H=M_PI/(csin(M_PI*U)*H);
      return F*H;
    } 
}


int fct(int f)
{
  if ( f < 0 ) 
    return 0;
  else if ( f == 0 ) 
    return 1;
  else
    return(f * fct(f - 1));
}


double complex wpsipg(double complex Z, int K)                                                  {    
  
  double complex U,V,H,R,P;                                      
  
  double C[6][5];
  static double DELTA = 5e-13;                                                
  static double R1 = 1.;
  static double HF = 1./2.;
  double C1 = pow(M_PI,2);
  double C2 = 2*pow(M_PI,3);
  double C3 = 2*pow(M_PI,4);
  double C4 = 8*pow(M_PI,5);
 
  C[0][0] =  8.33333333333333333e-2;
  C[1][0] = -8.33333333333333333e-3;
  C[2][0] =  3.96825396825396825e-3;
  C[3][0] = -4.16666666666666667e-3;
  C[4][0] =  7.57575757575757576e-3;
  C[5][0] = -2.10927960927960928e-2;
  
  C[0][1] =  1.66666666666666667e-1;
  C[1][1] = -3.33333333333333333e-2;
  C[2][1] =  2.38095238095238095e-2;
  C[3][1] = -3.33333333333333333e-2;
  C[4][1] =  7.57575757575757576e-2;
  C[5][1] = -2.53113553113553114e-1;
  
  C[0][2] =  5.00000000000000000e-1;
  C[1][2] = -1.66666666666666667e-1;
  C[2][2] =  1.66666666666666667e-1;
  C[3][2] = -3.00000000000000000e-1;
  C[4][2] =  8.33333333333333333e-1;
  C[5][2] = -3.29047619047619048e+0;
  
  C[0][3] =  2.00000000000000000e+0;
  C[1][3] = -1.00000000000000000e+0;
  C[2][3] =  1.33333333333333333e+0;
  C[3][3] = -3.00000000000000000e+0;
  C[4][3] =  1.00000000000000000e+1;
  C[5][3] = -4.60666666666666667e+1;

  C[0][4] =  1.00000000000000000e+1;
  C[1][4] = -7.00000000000000000e+0;
  C[2][4] =  1.20000000000000000e+1;
  C[3][4] = -3.30000000000000000e+1;
  C[4][4] =  1.30000000000000000e+2;
  C[5][4] = -6.91000000000000000e+2;
  
                                                                                
  U=Z;                                                                       
  double X = creal(U);                                                                       
  double A = fabs(X);                                                                  
  if(K < 0 || K > 4)                                       
    {       
      H=0;
      printf("Error");
      exit(1);      
    }
  else if(fabs(cimag(U)) < DELTA && fabs(X+nearbyint(A)) < DELTA)
    {
      H=0;
      printf("Error");
      exit(1);
    }
  else
    {                                                                      
      int K1=K+1;                                                                   
      if(X < 0) U=-U;
      V=U;                                                                
      H=0;
      if(A < 15)
        {                                                       
          H = cpow(V,-K1);
          int i;
          for(i = 1; i <= 14 - (int)A; i++)                                                    
            {
              V=V+1;
              H=H + cpow(V,-K1);
            }
          V=V+1;
        }
      R = cpow(V,-2);
      P = R*C[5][K];
      int i;
      for(i = 5; i >=1; i--)                                                          
        P=R*(C[i-1][K]+P);
                                                           
      H = -pow(-1,K)*(fct(K)*H+(V*(fct(K-1)+P)+HF*fct(K))*cpow(V,-K1));
      if(K == 0) H = H + clog(V);                                                  
      if(X < 0) 
        {                                                        
          V=M_PI*U;
          X=V;
          double Y=cimag(V);
          double A=sin(X);                                                         
          double B=cos(X);
          double T=tanh(Y);
          P=(B - I*A*T)/(A + I*B*T);
          if(K == 0) 
            {
              H=H+1/U+M_PI*P;
            }
          else if(K == 1)
            {
              H=-H + cpow(U,-2)+C1*(cpow(P,2)+1);
            }
          else if(K == 2)
            {
              H=H+2*cpow(U,-3) + C2*P*(cpow(P,2)+1);
            }
          else if(K == 3)
            {
              R=P*P;                                                               
              H=-H+6*cpow(U,-4) + C3*((3*R+4)*R+1);
            }
          else if(K == 4)
            {                                                   
              R=P*P;
              H=H+24*cpow(U,-5)+C4*P*((3*R+5)*R+2);
            }                           
        }                                                           
    }                                                               
  return H;                                                                   
}
              
#endif  /* __CGAMMA_H__ */
