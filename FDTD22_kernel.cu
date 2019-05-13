//--------------------------
// FDTD22
//--------------------------

#include <math.h>
#include <omp.h>

__global__ void FDTD22(
    float *Ex, float *Ey, float *Ez,
    float *Hx, float *Hy, float *Hz,
    float *CEx, float *CEy, float *CEz,
    float *CEx_dxyz_A, float *CEy_dxyz_A, float *CEz_dxyz_A,
    float *CH_dxyz_A,
    unsigned char *val,
    unsigned long Nx, unsigned long Ny, unsigned long Nz, unsigned long Nxy,
    double Ds, double dt
    )
{
  unsigned long X, Y, Z, ID;

  X = threadIdx.x + blockIdx.x * blockDim.x;
  Y = threadIdx.y + blockIdx.y * blockDim.y;
  Z = threadIdx.z + blockIdx.z * blockDim.z;
  ID = X + Y * Nx + Z * Nx * Ny;

  //--------------------------
  //解析空間の磁界計算
  //--------------------------

  // Magnetic Field Hx //
  if(          X < Nx     &&
      Y > 0 && Y < Ny - 2 &&
      Z > 0 && Z < Nz - 2)
  {
    Hx[ID] = Hx[ID]
      - CH_dxyz_A[val[ID]] * (Ez[ID+Nx] - Ez[ID])
      + CH_dxyz_A[val[ID]] * (Ey[ID+Nxy] - Ey[ID]);
  }
  // Magnetic Field Hy //
  if(X > 0 && X<Nx-2 &&
      Y<Ny   &&
      Z > 0 && Z<Nz-2)
  {
    Hy[ID] = Hy[ID]
      - CH_dxyz_A[val[ID]] * (Ex[ID+Nxy] - Ex[ID])
      + CH_dxyz_A[val[ID]] * (Ez[ID+1] - Ez[ID]);
  }
  // Magnetic Field Hz //
  if(X > 0 && X<Nx-2 &&
      Y > 0 && Y<Ny-2 &&
      Z<Nz)
  {
    Hz[ID] = Hz[ID]
      - CH_dxyz_A[val[ID]] * (Ey[ID+1] - Ey[ID])
      + CH_dxyz_A[val[ID]] * (Ex[ID+Nx] - Ex[ID]);
  }


  // Electric Field Ex //
  if(        X<Nx   &&
      Y>=2 && Y<Ny-2 &&
      Z>=2 && Z<Nz-2)
  {
    Ex[ID] = CEx[val[ID]] * Ex[ID]
      + CEx_dxyz_A[val[ID]] * (Hz[ID]   - Hz[ID-Nx])
      - CEx_dxyz_A[val[ID]] * (Hy[ID]   - Hy[ID-Nxy]);
  }
  // Electric Field Ey //
  if(X>=2 && X<Nx-2 &&
      Y<Ny   &&
      Z>=2 && Z<Nz-2)
  {
    Ey[ID] = CEy[val[ID]] * Ey[ID]
      + CEy_dxyz_A[val[ID]] * (Hx[ID]   - Hx[ID-Nxy])
      - CEy_dxyz_A[val[ID]] * (Hz[ID]   - Hz[ID-1]);
  }
  // Electric Field Ez //
  if(X>=2 && X<Nx-2 &&
      Y>=2 && Y<Ny-2 &&
      Z<Nz  )
  {
    Ez[ID] = CEz[val[ID]] * Ez[ID]
      + CEz_dxyz_A[val[ID]] * (Hy[ID]   - Hy[ID-1])
      - CEz_dxyz_A[val[ID]] * (Hx[ID]   - Hx[ID-Nx]);
  }
}

