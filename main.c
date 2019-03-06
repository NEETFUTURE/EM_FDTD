
#ifdef _MSC_VER//不要な警告表示の中止処理の追加
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "FDTD24_kernel.c"
#include "FDTD22_kernel.c"

// FDTD22法の場合はSCHEMEに2, 24法の場合は4を選択する
#if SCHEME==2 || SCHEME==4
#else
#error scheme symbol invalid
#endif

//Cell Size
#define Dx 0.04
#define Dy Dx
#define Dz Dx
#define Ds Dx

int main
(
 int    argc, // Argument Count
 char **argv  // Argument Vector
 )
{
  FILE *fp;
  char scheme[][10] = {"FDTD22", "FDTD24"};
  unsigned long i,j,k;

  unsigned long Nx = 101;
  unsigned long Ny = 101;
  unsigned long Nz = 102;
  unsigned long Nt = 2000;

  if(argc==5){
    Nx = (int)atoi(argv[1]);
    Ny = (int)atoi(argv[2]);
    Nz = (int)atoi(argv[3]);
    Nt = (int)atoi(argv[4]);
  }

  float *Ex;
  float *Ey;
  float *Ez;
  float *Ex_;
  float *Ey_;
  float *Ez_;

  float *Hx;
  float *Hy;
  float *Hz;

  float time_ini = 5.0e-9f; // Peak Time of Gaussian Pulse

  unsigned long Nxy = Nx*Ny;

  unsigned long elements = Nx*Ny*Nz*sizeof(float);

  printf("Scheme is %s\n", scheme[SCHEME-1]);
  printf("Nx = %lu\n", Nx);
  printf("Ny = %lu\n", Ny);
  printf("Nz = %lu\n", Nz);
  printf("elements = %lu\n", elements);

  // 媒質係数のインデックス
  unsigned char *val;

  float *CEx;
  float *CEy;
  float *CEz;
  float *CEx_dxyz_A;
  float *CEx_dxyz_B;
  float *CEy_dxyz_A;
  float *CEy_dxyz_B;
  float *CEz_dxyz_A;
  float *CEz_dxyz_B;
  float *CH_dxyz_A;
  float *CH_dxyz_B;


  double CG1 = 1.0e+0;
  double CG2 = 1.5e+18;

#ifdef DEBUG
  printf("DEBUG MODE!\n");
#endif

  printf("start malloc!\n");

  Ex = (float*)malloc(elements);
  Ey = (float*)malloc(elements);
  Ez = (float*)malloc(elements);
  Hx = (float*)malloc(elements);
  Hy = (float*)malloc(elements);
  Hz = (float*)malloc(elements);

  val = (unsigned char *)malloc(Nx*Ny*Nz*sizeof(unsigned char));

  CEx = (float*)malloc(2*sizeof(float));
  CEy = (float*)malloc(2*sizeof(float));
  CEz = (float*)malloc(2*sizeof(float));

  CEx_dxyz_A = (float*)malloc(2*sizeof(float));
  CEx_dxyz_B = (float*)malloc(2*sizeof(float));
  CEy_dxyz_A = (float*)malloc(2*sizeof(float));
  CEy_dxyz_B = (float*)malloc(2*sizeof(float));
  CEz_dxyz_A = (float*)malloc(2*sizeof(float));
  CEz_dxyz_B = (float*)malloc(2*sizeof(float));
  CH_dxyz_A  = (float*)malloc(2*sizeof(float));
  CH_dxyz_B  = (float*)malloc(2*sizeof(float));

  printf("finish malloc!\n");

  //-----------------
  //媒質パラメータ設定
  //-----------------

  double e[3]; //誘電率（ε）
  double u[3]; //透磁率（μ）
  double s[3]; //導電率（σ)
  double pi = acos(-1.0);

  //真空
  e[0] = 8.85418782e-12;
  u[0] = pi * 4.e-7;
  s[0] = 0.0;
  e[1] = 8.85418782e-12;
  //コンクリート
  e[1]  = 15.0 * e[0];
  u[1]  = u[0];
  s[1]  = 0.015;

  double CFL = 0.13105f;//クーラン数(0 < CFL <= 0.7)
  double c0 = 1.0 / sqrt(e[0] * u[0]);// = 2.99792458e+8f;//光速[m/s]
  double InitialTimeStep = 1.0 /c0 * Dx;
  double dt = CFL * InitialTimeStep;
  printf("dt = %g\n", dt);

  for(i=0;i<Nx;i++)
  {
    for(j=0;j<Ny;j++)
    {
      for(k=0;k<Nz;k++)
      {
        unsigned long ID = i+j*Nx+k*Nxy;

        val[ID] = 0;

        if( (i>=0 && i<=9) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          val[ID] = 1;
        }
        if( (i>=Nx-1-10 && i<=Nx-1) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=0 && j<=9) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=Ny-1-10 && j<=Ny-1) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=0 && k<=9) )
        {
          val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=Nz-1-10 && k<=Nz-1) )
        {
          val[ID] = 1;
        }
      }
    }
  }

  for(i=0;i<2;i++)
  {
    CEx[i]   = (1.0 - (s[i]*dt)/(2.0*e[i])) / (1.0 + (s[i]*dt)/(2.0*e[i]));
    CEy[i]   = (1.0 - (s[i]*dt)/(2.0*e[i])) / (1.0 + (s[i]*dt)/(2.0*e[i]));
    CEz[i]   = (1.0 - (s[i]*dt)/(2.0*e[i])) / (1.0 + (s[i]*dt)/(2.0*e[i]));


#if SCHEME == 2
    CEx_dxyz_A[i] = (dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CEy_dxyz_A[i] = (dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CEz_dxyz_A[i] = (dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CH_dxyz_A[i]  = dt / u[i] / Ds;
#elif SCHEME == 4
    CEx_dxyz_A[i] = 9.0/8.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CEy_dxyz_A[i] = 9.0/8.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CEz_dxyz_A[i] = 9.0/8.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CEx_dxyz_B[i] = 1.0/24.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CEy_dxyz_B[i] = 1.0/24.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CEz_dxyz_B[i] = 1.0/24.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    CH_dxyz_A[i]  =  9.0 /  8.0 * dt / u[i] / Ds;
    CH_dxyz_B[i]  =  1.0 / 24.0 * dt / u[i] / Ds;
#endif
  }


  for(unsigned long ID=0; ID<Nx*Ny*Nz; ID++)
  {
    Ex[ID] = 0.0;
    Ey[ID] = 0.0;
    Ez[ID] = 0.0;
    Hx[ID] = 0.0;
    Hy[ID] = 0.0;
    Hz[ID] = 0.0;
  }


  unsigned long ant_pos_x = Nx/2;
  unsigned long ant_pos_y = Ny/2;
  unsigned long ant_pos_z = Nz/2;

  //アンテナの座標
  printf("ant_pos_x = %lu\n", ant_pos_x);
  printf("ant_pos_y = %lu\n", ant_pos_y);
  printf("ant_pos_z = %lu\n", ant_pos_z);



#if SCHEME==2
  FDTD22(
      Ex,  Ey,  Ez,
      Hx,  Hy,  Hz,
      CEx, CEy, CEz,
      CEx_dxyz_A, CEy_dxyz_A, CEz_dxyz_A,
      CH_dxyz_A,
      val,
      Nx, Ny, Nz, Nt,
      Ds, dt,
      CG1, CG2,
      ant_pos_x,
      ant_pos_y,
      ant_pos_z,
      time_ini
      );
#elif SCHEME==1
  FDTD24(
      Ex,  Ey,  Ez,
      Hx,  Hy,  Hz,
      CEx, CEy, CEz,
      CEx_dxyz_A, CEy_dxyz_A, CEz_dxyz_A,
      CEx_dxyz_B, CEy_dxyz_B, CEz_dxyz_B,
      CH_dxyz_A,
      CH_dxyz_B,
      val,
      Nx, Ny, Nz, Nt,
      Ds, dt,
      CG1, CG2,
      ant_pos_x,
      ant_pos_y,
      ant_pos_z,
      time_ini
      );
#endif

  time_t now = time(NULL);
  struct tm *pnow = localtime(&now);
  char filename[80];

  // Z方向の電場をファイル出力
  sprintf(filename, "%04d%02d%02d%02d%02d%02d_EM_%s_Ez.bin",pnow->tm_year+1900,pnow->tm_mon+1,pnow->tm_mday,pnow->tm_hour,pnow->tm_min,pnow->tm_sec, scheme[(int)(SCHEME/2)-1]);
  fp = fopen(filename,"wb");
  fwrite(&Ez[Nxy*Nz/2], sizeof(float), Nx*Ny, fp);
  fclose(fp);

  return 0;
}
