//--------------------------
// FDTD24
//--------------------------

#include <math.h>

void FDTD24(
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
    )
{
  unsigned long i,j,k;
  unsigned long Nxy = Nx*Ny;
  float wave;

  for(int step=0; step<Nt; step++)
  {

#ifdef DEBUG
    if(step%10==0)
    {
      printf("step = %d\n", step);
    }
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
#pragma omp parallel for private(i,j)
    for(k=0; k<Nz; k++)
    {
      for(j=0; j<Ny; j++)
      {
        for(i=0; i<Nx; i++)
        {

          unsigned long ID = i+j*Nx+k*Nxy;

          // Magnetic Field Hx //
          if(         i<Nx   &&
              j > 0 && j<Ny-2 &&
              k > 0 && k<Nz-2)
          {
            Hx[ID] = Hx[ID]
              - CH_dxyz_A[val[ID]] * (Ez[ID+Nx] - Ez[ID])
              + CH_dxyz_B[val[ID]] * (Ez[ID+2*Nx] - Ez[ID-Nx])
              + CH_dxyz_A[val[ID]] * (Ey[ID+Nxy] - Ey[ID])
              - CH_dxyz_B[val[ID]] * (Ey[ID+2*Nxy] - Ey[ID-Nxy]);
          }
          // Magnetic Field Hy //
          if(i > 0 && i<Nx-2 &&
              j<Ny   &&
              k > 0 && k<Nz-2)
          {
            Hy[ID] = Hy[ID]
              - CH_dxyz_A[val[ID]] * (Ex[ID+Nxy] - Ex[ID])
              + CH_dxyz_B[val[ID]] * (Ex[ID+2*Nxy] - Ex[ID-Nxy])
              + CH_dxyz_A[val[ID]] * (Ez[ID+1] - Ez[ID])
              - CH_dxyz_B[val[ID]] * (Ez[ID+2] - Ez[ID-1]);
          }
          // Magnetic Field Hz //
          if(i > 0 && i<Nx-2 &&
              j > 0 && j<Ny-2 &&
              k<Nz)
          {
            Hz[ID] = Hz[ID]
              - CH_dxyz_A[val[ID]] * (Ey[ID+1] - Ey[ID])
              + CH_dxyz_B[val[ID]] * (Ey[ID+2] - Ey[ID-1])
              + CH_dxyz_A[val[ID]] * (Ex[ID+Nx] - Ex[ID])
              - CH_dxyz_B[val[ID]] * (Ex[ID+2*Nx] - Ex[ID-Nx]);
          }
        }
      }
    }

    //--------------------------
    //解析空間の電界計算
    //--------------------------
#pragma omp parallel for private(i,j)
    for(k = 0; k < Nz; k++)
    {
      for(j = 0; j < Ny; j++)
      {
        for(i = 0; i < Nx; i++)
        {

          unsigned long ID = i+j*Nx+k*Nxy;

          // Electric Field Ex //
          if(        i<Nx   &&
              j>=2 && j<Ny-2 &&
              k>=2 && k<Nz-2)
          {
            Ex[ID] = CEx[val[ID]] * Ex[ID]
              + CEx_dxyz_A[val[ID]] * (Hz[ID]   - Hz[ID-Nx])
              - CEx_dxyz_B[val[ID]] * (Hz[ID+Nx] - Hz[ID-2*Nx])
              - CEx_dxyz_A[val[ID]] * (Hy[ID]   - Hy[ID-Nxy])
              + CEx_dxyz_B[val[ID]] * (Hy[ID+Nxy] - Hy[ID-2*Nxy]);
          }
          // Electric Field Ey //
          if(i>=2 && i<Nx-2 &&
              j<Ny   &&
              k>=2 && k<Nz-2)
          {
            Ey[ID] = CEy[val[ID]] * Ey[ID]
              + CEy_dxyz_A[val[ID]] * (Hx[ID]   - Hx[ID-Nxy])
              - CEy_dxyz_B[val[ID]] * (Hx[ID+Nxy] - Hx[ID-2*Nxy])
              - CEy_dxyz_A[val[ID]] * (Hz[ID]   - Hz[ID-1])
              + CEy_dxyz_B[val[ID]] * (Hz[ID+1] - Hz[ID-2]);
          }
          // Electric Field Ez //
          if(i>=2 && i<Nx-2 &&
              j>=2 && j<Ny-2 &&
              k<Nz  )
          {
            Ez[ID] = CEz[val[ID]] * Ez[ID]
              + CEz_dxyz_A[val[ID]] * (Hy[ID]   - Hy[ID-1])
              - CEz_dxyz_B[val[ID]] * (Hy[ID+1] - Hy[ID-2])
              - CEz_dxyz_A[val[ID]] * (Hx[ID]   - Hx[ID-Nx])
              + CEz_dxyz_B[val[ID]] * (Hx[ID+Nx] - Hx[ID-2*Nx]);
          }
        }
      }
    }
  }
}

