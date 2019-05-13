
#ifdef _MSC_VER//不要な警告表示の中止処理の追加
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <timer.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#define HTH cudaMemcpyHostToHost
#define HTD cudaMemcpyHostToDevice
#define DTH cudaMemcpyDeviceToHost
#define DTD cudaMemcpyDeviceToDevice

#define Db_x 32
#define Db_y 1
#define Db_z 8

#define Dg_x (Nx / Db_x)
#define Dg_y (Ny / Db_y)
#define Dg_z (Nz / Db_z)

#define THRESHOLD 200

#include "Antenna.cu"
#if SCHEME==2
#include "FDTD22_kernel.cu"
#elif SCHEME==4
#include "FDTD24_kernel.cu"
#endif

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
  unsigned long YBF = 8;

  if(argc>=4){
    Nx = (int)atoi(argv[1]);
    Ny = (int)atoi(argv[2]);
    Nz = (int)atoi(argv[3]);
  }
  if(argc>=5){
    Nt = (int)atoi(argv[4]);
  }
  if(argc>=6){
    YBF = (unsigned long)atoi(argv[5]);
  }

  float *Ex;
  float *Ey;
  float *Ez;
  float *H_Ez;
  float *Ex_;
  float *Ey_;
  float *Ez_;

  float *Hx;
  float *Hy;
  float *Hz;

  float time_ini = 5.0e-9f; // Peak Time of Gaussian Pulse

  unsigned long Nxy = Nx*Ny;

  unsigned long elements = Nx*Ny*Nz*sizeof(float);

  dim3 Db_FDTD(Db_x, Db_y, Db_z);
  dim3 Dg_FDTD(Dg_x, Dg_y, Dg_z);
  dim3 Db_Antenna(   16,    16,  1);
  dim3 Dg_Antenna(Nx/16, Ny/16, Nz);

  printf("Scheme is %s\n", scheme[(int)(SCHEME/2-1)]);
  printf("Nx = %lu\n", Nx);
  printf("Ny = %lu\n", Ny);
  printf("Nz = %lu\n", Nz);
  printf("elements = %lu\n", elements);
	printf("Thread Block: (%d, %d, %d)\n", Db_x, Db_y, Db_z);
	printf("Grid Block  : (%d, %d, %d)\n", Dg_x, Dg_y, Dg_z);

  // デバイス媒質係数のインデックス
  unsigned char *val;

  // ホスト媒質係数のインデックス
  unsigned char *H_val;

  // デバイス側の媒質係数変数
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

  // ホスト側の媒質係数変数
  float *H_CEx;
  float *H_CEy;
  float *H_CEz;
  float *H_CEx_dxyz_A;
  float *H_CEx_dxyz_B;
  float *H_CEy_dxyz_A;
  float *H_CEy_dxyz_B;
  float *H_CEz_dxyz_A;
  float *H_CEz_dxyz_B;
  float *H_CH_dxyz_A;
  float *H_CH_dxyz_B;

  float wave;
  double CG1 = 1.0e+0;
  double CG2 = 1.5e+18;

  double processing_time;
  StopWatchInterface *timer = NULL;

#ifdef DEBUG
  printf("DEBUG MODE!\n");
#endif

  printf("start malloc!\n");

  // デバイス側配列の確保
  checkCudaErrors(cudaMalloc((void **)&Ex, elements));
  checkCudaErrors(cudaMalloc((void **)&Ey, elements));
  checkCudaErrors(cudaMalloc((void **)&Ez, elements));
  checkCudaErrors(cudaMalloc((void **)&Hx, elements));
  checkCudaErrors(cudaMalloc((void **)&Hy, elements));
  checkCudaErrors(cudaMalloc((void **)&Hz, elements));

  checkCudaErrors(cudaMalloc((void **)&val, sizeof(unsigned char)*Nx*Ny*Nz));

  checkCudaErrors(cudaMalloc((void **)&CEx, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CEy, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CEz, sizeof(float)*2));

  checkCudaErrors(cudaMalloc((void **)&CEx_dxyz_A, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CEy_dxyz_A, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CEz_dxyz_A, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CH_dxyz_A, sizeof(float)*2));
#if SCHEME == 4
  checkCudaErrors(cudaMalloc((void **)&CEx_dxyz_B, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CEy_dxyz_B, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CEz_dxyz_B, sizeof(float)*2));
  checkCudaErrors(cudaMalloc((void **)&CH_dxyz_B, sizeof(float)*2));
#endif

  // デバイス側配列の初期化
  checkCudaErrors(cudaMemset(Ex, 0, elements));
  checkCudaErrors(cudaMemset(Ey, 0, elements));
  checkCudaErrors(cudaMemset(Ez, 0, elements));
  checkCudaErrors(cudaMemset(Hx, 0, elements));
  checkCudaErrors(cudaMemset(Hy, 0, elements));
  checkCudaErrors(cudaMemset(Hz, 0, elements));

  checkCudaErrors(cudaMemset(val, 0, sizeof(unsigned char)*Nx*Ny*Nz));

  checkCudaErrors(cudaMemset(CEx, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CEy, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CEz, 0, sizeof(float)*2));

  checkCudaErrors(cudaMemset(CEx_dxyz_A, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CEy_dxyz_A, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CEz_dxyz_A, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CH_dxyz_A, 0, sizeof(float)*2));
#if SCHEME == 4
  checkCudaErrors(cudaMemset(CEx_dxyz_B, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CEy_dxyz_B, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CEz_dxyz_B, 0, sizeof(float)*2));
  checkCudaErrors(cudaMemset(CH_dxyz_B, 0, sizeof(float)*2));
#endif

  // ホスト側配列の確保
  checkCudaErrors(cudaHostAlloc((void **)&H_Ez, elements, cudaHostAllocPortable));

  H_val = (unsigned char *)malloc(Nx*Ny*Nz*sizeof(unsigned char));

  H_CEx = (float*)malloc(2*sizeof(float));
  H_CEy = (float*)malloc(2*sizeof(float));
  H_CEz = (float*)malloc(2*sizeof(float));

  H_CEx_dxyz_A = (float*)malloc(2*sizeof(float));
  H_CEx_dxyz_B = (float*)malloc(2*sizeof(float));
  H_CEy_dxyz_A = (float*)malloc(2*sizeof(float));
  H_CEy_dxyz_B = (float*)malloc(2*sizeof(float));
  H_CEz_dxyz_A = (float*)malloc(2*sizeof(float));
  H_CEz_dxyz_B = (float*)malloc(2*sizeof(float));
  H_CH_dxyz_A  = (float*)malloc(2*sizeof(float));
  H_CH_dxyz_B  = (float*)malloc(2*sizeof(float));

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

        H_val[ID] = 0;

        if( (i>=0 && i<=Nx/6) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          H_val[ID] = 1;
        }
        if( (i>=Nx-1-Nx/6 && i<=Nx-1) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          H_val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=0 && j<=Ny/6) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          H_val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=Ny-1-Ny/6 && j<=Ny-1) && (k>=Nz/4 && k<=Nz*3/4) )
        {
          H_val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=0 && k<=Nz/6) )
        {
          H_val[ID] = 1;
        }
        if( (i>=Nx/4 && i<=Nx*3/4-1) && (j>=Ny/4 && j<=Ny*3/4-1) && (k>=Nz-1-Nz/6 && k<=Nz-1) )
        {
          H_val[ID] = 1;
        }
      }
    }
  }
  checkCudaErrors(cudaMemcpy(val, H_val, sizeof(unsigned char)*Nx*Ny*Nz, HTD));

  for(i=0;i<2;i++)
  {
    H_CEx[i]   = (1.0 - (s[i]*dt)/(2.0*e[i])) / (1.0 + (s[i]*dt)/(2.0*e[i]));
    H_CEy[i]   = (1.0 - (s[i]*dt)/(2.0*e[i])) / (1.0 + (s[i]*dt)/(2.0*e[i]));
    H_CEz[i]   = (1.0 - (s[i]*dt)/(2.0*e[i])) / (1.0 + (s[i]*dt)/(2.0*e[i]));
    checkCudaErrors(cudaMemcpy(CEx, H_CEx, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEy, H_CEy, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEz, H_CEz, sizeof(float)*2, HTD));

#if SCHEME == 2
    H_CEx_dxyz_A[i] = (dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CEy_dxyz_A[i] = (dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CEz_dxyz_A[i] = (dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CH_dxyz_A[i]  = dt / u[i] / Ds;
    checkCudaErrors(cudaMemcpy(CEx_dxyz_A, H_CEx_dxyz_A, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEy_dxyz_A, H_CEy_dxyz_A, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEz_dxyz_A, H_CEz_dxyz_A, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CH_dxyz_A, H_CH_dxyz_A, sizeof(float)*2, HTD));
#elif SCHEME == 4
    H_CEx_dxyz_A[i] = 9.0/8.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CEy_dxyz_A[i] = 9.0/8.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CEz_dxyz_A[i] = 9.0/8.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CEx_dxyz_B[i] = 1.0/24.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CEy_dxyz_B[i] = 1.0/24.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CEz_dxyz_B[i] = 1.0/24.0*(dt/e[i])/(1.0+(s[i]*dt)/(2.0*e[i]))/Ds;
    H_CH_dxyz_A[i]  =  9.0 /  8.0 * dt / u[i] / Ds;
    H_CH_dxyz_B[i]  =  1.0 / 24.0 * dt / u[i] / Ds;
    checkCudaErrors(cudaMemcpy(CEx_dxyz_A, H_CEx_dxyz_A, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEy_dxyz_A, H_CEy_dxyz_A, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEz_dxyz_A, H_CEz_dxyz_A, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEx_dxyz_B, H_CEx_dxyz_B, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEy_dxyz_B, H_CEy_dxyz_B, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CEz_dxyz_B, H_CEz_dxyz_B, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CH_dxyz_A, H_CH_dxyz_A, sizeof(float)*2, HTD));
    checkCudaErrors(cudaMemcpy(CH_dxyz_B, H_CH_dxyz_B, sizeof(float)*2, HTD));
#endif
  }


  unsigned long ant_pos_x = Nx/2;
  unsigned long ant_pos_y = Ny/2;
  unsigned long ant_pos_z = Nz/2;

  //アンテナの座標
  printf("ant_pos_x = %lu\n", ant_pos_x);
  printf("ant_pos_y = %lu\n", ant_pos_y);
  printf("ant_pos_z = %lu\n", ant_pos_z);


  sdkCreateTimer(&timer);
  sdkResetTimer(&timer);

  printf("\nCalculation Start\n");


  sdkStartTimer(&timer);

  for(int step=0; step < Nt+1; step++)
  {
#ifdef DEBUG
    if(step % 100 == 0)
      printf("step = %d\n", step);
#endif

    wave = (float)(CG1 * expf(-((((float)step * dt) - time_ini) * (((float)step * dt) - time_ini)) * CG2))/Ds;


    Antenna <<< Dg_Antenna, Db_Antenna >>> (
        Ez, Nx, Ny, Nz, wave,
        ant_pos_x, ant_pos_y, ant_pos_z);

#if SCHEME==2
    FDTD22 <<< Dg_FDTD, Db_FDTD >>> (
        Ex,  Ey,  Ez,
        Hx,  Hy,  Hz,
        CEx, CEy, CEz,
        CEx_dxyz_A, CEy_dxyz_A, CEz_dxyz_A,
        CH_dxyz_A,
        val,
        Nx, Ny, Nz, Nxy,
        Ds, dt
        );
#elif SCHEME==4
    FDTD24 <<< Dg_FDTD, Db_FDTD >>> (
        Ex,  Ey,  Ez,
        Hx,  Hy,  Hz,
        CEx, CEy, CEz,
        CEx_dxyz_A, CEy_dxyz_A, CEz_dxyz_A,
        CEx_dxyz_B, CEy_dxyz_B, CEz_dxyz_B,
        CH_dxyz_A,
        CH_dxyz_B,
        val,
        Nx, Ny, Nz, Nt,
        Ds, dt
        );
#endif
    checkCudaErrors(cudaThreadSynchronize());

  }

  sdkStopTimer(&timer);
  processing_time = sdkGetTimerValue(&timer);
  printf("\nCalculation End\n");
	printf("\nProcessing Time : %f [ms]\n", processing_time);

  char cfg_str[300];
  char ybf_str[80];

  sprintf(cfg_str, "_CUDA");

  time_t now = time(NULL);
  struct tm *pnow = localtime(&now);
  char filename[80];

#ifdef DEBUG
  checkCudaErrors(cudaMemcpy(H_Ez, Ez, elements, DTH));
  // Z方向の電場をファイル出力
  sprintf(filename, "%04d%02d%02d%02d%02d%02d_EM_%s%s_Ez.bin",pnow->tm_year+1900,pnow->tm_mon+1,pnow->tm_mday,pnow->tm_hour,pnow->tm_min,pnow->tm_sec, scheme[(int)(SCHEME/2)-1], cfg_str);
  fp = fopen(filename,"wb");
  fwrite(&H_Ez[Nxy*Nz/2], sizeof(float), Nx*Ny, fp);
  fclose(fp);
#endif

  // 計算時間の記録
#ifdef DEBUG_LAP
  sprintf(filename, "EM_%s_%lu_%lu%s_timelap.csv", scheme[(int)(SCHEME/2)-1], Nx, Nt, cfg_str);
  fp = fopen(filename,"w");
  for(int step=0; step < Nt; step++)
  {
    fprintf(fp, "%lf\n", laptime[step]);
  }
  fclose(fp);
#else
  sprintf(filename, "EM_%s_%lu_%lu%s_Ez.csv", scheme[(int)(SCHEME/2)-1], Nx, Nt-THRESHOLD, cfg_str);
  fp = fopen(filename,"a");
  fprintf(fp, "%lf\n", processing_time);
  fclose(fp);
#endif


  sdkDeleteTimer(&timer);
  return 0;
}
