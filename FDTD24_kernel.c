//--------------------------
// FDTD24
//--------------------------

#include <math.h>
#include <omp.h>

#ifdef DEBUG_LAP
double* FDTD24(
#else
double FDTD24(
#endif
    float *Ex, float *Ey, float *Ez,
    float *Hx, float *Hy, float *Hz,
    float *CEx, float *CEy, float *CEz,
    float *CEx_dxyz_A, float *CEy_dxyz_A, float *CEz_dxyz_A,
    float *CEx_dxyz_B, float *CEy_dxyz_B, float *CEz_dxyz_B,
    float *CH_dxyz_A,
    float *CH_dxyz_B,
    unsigned char *val,
    unsigned long Nx, unsigned long Ny, unsigned long Nz, int Nt, double Ds, double dt,
    double CG1, double CG2,
    unsigned long ant_pos_x, unsigned long ant_pos_y, unsigned long ant_pos_z, float time_ini
#ifdef TILE
    ,unsigned long YBF
#endif
    )
{
  double start,end;
  unsigned long i,j,k;
  unsigned long Nxy = Nx*Ny;
  float wave;
#ifdef DEBUG_LAP
  double lap_start = 0.0;
  double *laptimes;
  laptimes = (double *)malloc(sizeof(double)*Nt);
#endif

  start = omp_get_wtime();
  for(int step=0; step<Nt+1; step++)
  {
    if(step == THRESHOLD)
    {
      start = omp_get_wtime();
    }

#ifdef DEBUG
    if(step%10==0)
    {
      printf("step = %d\n", step);
    }
#endif
#ifdef DEBUG_LAP
    lap_start = omp_get_wtime();
#endif

    wave = (float)(CG1 * expf(-((((float)step * dt) - time_ini) * (((float)step * dt) - time_ini)) * CG2));

    for(unsigned long k=ant_pos_z-2; k<=ant_pos_z+2; k++)
    {
      Ez[ant_pos_x+ant_pos_y*Nx+k*Nxy] = 0.0;
    }
    Ez[ant_pos_x+ant_pos_y*Nx+ant_pos_z*Nxy] = wave/Ds;

    //--------------------------
    //解析空間の磁界計算
    //--------------------------

    /*                   */
    /* Magnetic Field Hx */
    /*                   */
#ifdef COLLAPSE
#pragma omp parallel for collapse(2)
#else
#pragma omp parallel for
#endif

#ifdef TILE
    for(unsigned long jj=1; jj < Ny-2; jj += YBF)
    {
#endif
    for(unsigned long k=1; k<Nz-2; k++)
    {
#ifdef TILE
      unsigned long jmax = jj + YBF;
      if(jmax > Ny-2) jmax = Ny-2;
      for(unsigned long j = jj; j < jmax; j++)
#else
      for(unsigned long j=1; j<Ny-2; j++)
#endif
      {
#ifdef SIMD
#pragma omp simd
#endif
        for(unsigned long i=0; i<Nx; i++)
        {
          unsigned long ID = i + j * Nx + k * Nxy;

          // Magnetic Field Hx //
          //if(         i<Nx   &&
          //    j > 0 && j<Ny-2 &&
          //    k > 0 && k<Nz-2)
          //{
            Hx[ID] = Hx[ID]
              - CH_dxyz_A[val[ID]] * (Ez[ID+Nx] - Ez[ID])
              + CH_dxyz_B[val[ID]] * (Ez[ID+2*Nx] - Ez[ID-Nx])
              + CH_dxyz_A[val[ID]] * (Ey[ID+Nxy] - Ey[ID])
              - CH_dxyz_B[val[ID]] * (Ey[ID+2*Nxy] - Ey[ID-Nxy]);
          //}

        }
      }
    }
#ifdef TILE
    }
#endif


    /*                   */
    /* Magnetic Field Hy */
    /*                   */
#ifdef COLLAPSE
#pragma omp parallel for collapse(2)
#else
#pragma omp parallel for
#endif

#ifdef TILE
    for(unsigned long jj=0; jj < Ny; jj += YBF)
    {
#endif
    for(unsigned long k=1; k<Nz-2; k++)
    {
#ifdef TILE
      unsigned long jmax = jj + YBF;
      if(jmax > Ny-2) jmax = Ny-2;
      for(unsigned long j = jj; j < jmax; j++)
#else
      for(unsigned long j=0; j<Ny; j++)
#endif
      {
#ifdef SIMD
#pragma omp simd
#endif
        for(unsigned long i=1; i<Nx-2; i++)
        {
          unsigned long ID = i + j * Nx + k * Nxy;
          // Magnetic Field Hy //
          //if(i > 0 && i<Nx-2 &&
          //    j<Ny   &&
          //    k > 0 && k<Nz-2)
          //{
            Hy[ID] = Hy[ID]
              - CH_dxyz_A[val[ID]] * (Ex[ID+Nxy] - Ex[ID])
              + CH_dxyz_B[val[ID]] * (Ex[ID+2*Nxy] - Ex[ID-Nxy])
              + CH_dxyz_A[val[ID]] * (Ez[ID+1] - Ez[ID])
              - CH_dxyz_B[val[ID]] * (Ez[ID+2] - Ez[ID-1]);
          //}

        }
      }
    }
#ifdef TILE
    }
#endif


    /*                   */
    /* Magnetic Field Hz */
    /*                   */
#ifdef COLLAPSE
#pragma omp parallel for collapse(2)
#else
#pragma omp parallel for
#endif

#ifdef TILE
    for(unsigned long jj=1; jj < Ny-2; jj += YBF)
    {
#endif
    for(unsigned long k=0; k<Nz; k++)
    {
#ifdef TILE
      unsigned long jmax = jj + YBF;
      if(jmax > Ny-2) jmax = Ny-2;
      for(unsigned long j = jj; j < jmax; j++)
#else
      for(unsigned long j=1; j<Ny-2; j++)
#endif
      {
#ifdef SIMD
#pragma omp simd
#endif
        for(unsigned long i=1; i<Nx-2; i++)
        {
          unsigned long ID = i + j * Nx + k * Nxy;
          // Magnetic Field Hz //
          //if(i > 0 && i<Nx-2 &&
          //    j > 0 && j<Ny-2 &&
          //    k<Nz)
          //{
            Hz[ID] = Hz[ID]
              - CH_dxyz_A[val[ID]] * (Ey[ID+1] - Ey[ID])
              + CH_dxyz_B[val[ID]] * (Ey[ID+2] - Ey[ID-1])
              + CH_dxyz_A[val[ID]] * (Ex[ID+Nx] - Ex[ID])
              - CH_dxyz_B[val[ID]] * (Ex[ID+2*Nx] - Ex[ID-Nx]);
          //}

        }
      }
    }
#ifdef TILE
    }
#endif



    //--------------------------
    //解析空間の電界計算
    //--------------------------

    /*                   */
    /* Electric Field Ex */
    /*                   */
#ifdef COLLAPSE
#pragma omp parallel for collapse(2)
#else
#pragma omp parallel for
#endif

#ifdef TILE
    for(unsigned long jj=2; jj < Ny-2; jj += YBF)
    {
#endif
    for(unsigned long k = 2; k < Nz-2; k++)
    {
#ifdef TILE
      unsigned long jmax = jj + YBF;
      if(jmax > Ny-2) jmax = Ny-2;
      for(unsigned long j = jj; j < jmax; j++)
#else
      for(unsigned long j = 2; j < Ny-2; j++)
#endif
      {
#ifdef SIMD
#pragma omp simd
#endif
        for(unsigned long i = 0; i < Nx; i++)
        {
          unsigned long ID = i + j * Nx + k * Nxy;

          // Electric Field Ex //
          //if(        i<Nx   &&
          //    j>=2 && j<Ny-2 &&
          //    k>=2 && k<Nz-2)
          //{
            Ex[ID] = CEx[val[ID]] * Ex[ID]
              + CEx_dxyz_A[val[ID]] * (Hz[ID]   - Hz[ID-Nx])
              - CEx_dxyz_B[val[ID]] * (Hz[ID+Nx] - Hz[ID-2*Nx])
              - CEx_dxyz_A[val[ID]] * (Hy[ID]   - Hy[ID-Nxy])
              + CEx_dxyz_B[val[ID]] * (Hy[ID+Nxy] - Hy[ID-2*Nxy]);
          //}

        }
      }
    }
#ifdef TILE
    }
#endif


    /*                   */
    /* Electric Field Ey */
    /*                   */
#ifdef COLLAPSE
#pragma omp parallel for collapse(2)
#else
#pragma omp parallel for
#endif

#ifdef TILE
    for(unsigned long jj=0; jj < Ny; jj += YBF)
    {
#endif
    for(unsigned long k = 2; k < Nz-2; k++)
    {
#ifdef TILE
      unsigned long jmax = jj + YBF;
      if(jmax > Ny-2) jmax = Ny-2;
      for(unsigned long j = jj; j < jmax; j++)
#else
      for(unsigned long j = 0; j < Ny; j++)
#endif
      {
#ifdef SIMD
#pragma omp simd
#endif
        for(unsigned long i = 2; i < Nx-2; i++)
        {
          unsigned long ID = i + j * Nx + k * Nxy;

          // Electric Field Ey //
          //if(i>=2 && i<Nx-2 &&
          //    j<Ny   &&
          //    k>=2 && k<Nz-2)
          //{
            Ey[ID] = CEy[val[ID]] * Ey[ID]
              + CEy_dxyz_A[val[ID]] * (Hx[ID]   - Hx[ID-Nxy])
              - CEy_dxyz_B[val[ID]] * (Hx[ID+Nxy] - Hx[ID-2*Nxy])
              - CEy_dxyz_A[val[ID]] * (Hz[ID]   - Hz[ID-1])
              + CEy_dxyz_B[val[ID]] * (Hz[ID+1] - Hz[ID-2]);
          //}

        }
      }
    }
#ifdef TILE
    }
#endif


    /*                   */
    /* Electric Field Ez */
    /*                   */
#ifdef COLLAPSE
#pragma omp parallel for collapse(2)
#else
#pragma omp parallel for
#endif

#ifdef TILE
    for(unsigned long jj=2; jj < Ny-2; jj += YBF)
    {
#endif
    for(unsigned long k = 0; k < Nz; k++)
    {
#ifdef TILE
      unsigned long jmax = jj + YBF;
      if(jmax > Ny-2) jmax = Ny-2;
      for(unsigned long j = jj; j < jmax; j++)
#else
      for(unsigned long j = 2; j < Ny-2; j++)
#endif
      {
#ifdef SIMD
#pragma omp simd
#endif
        for(unsigned long i = 2; i < Nx-2; i++)
        {
          unsigned long ID = i + j * Nx + k * Nxy;

          // Electric Field Ez //
          //if(i>=2 && i<Nx-2 &&
          //    j>=2 && j<Ny-2 &&
          //    k<Nz  )
          //{
            Ez[ID] = CEz[val[ID]] * Ez[ID]
              + CEz_dxyz_A[val[ID]] * (Hy[ID]   - Hy[ID-1])
              - CEz_dxyz_B[val[ID]] * (Hy[ID+1] - Hy[ID-2])
              - CEz_dxyz_A[val[ID]] * (Hx[ID]   - Hx[ID-Nx])
              + CEz_dxyz_B[val[ID]] * (Hx[ID+Nx] - Hx[ID-2*Nx]);
          //}

        }
      }
    }
#ifdef TILE
    }
#endif


#ifdef DEBUG_LAP
    laptimes[step] = omp_get_wtime() - lap_start;
#endif
  }

#ifdef DEBUG_LAP
  return laptimes;
#else
  end = omp_get_wtime();
  return end - start;
#endif
}

