
__global__ void Antenna
(
 float *Ez,
 unsigned long Nx,
 unsigned long Ny,
 unsigned long Nz,
 float wave,
 unsigned long ant_pos_x, unsigned long ant_pos_y, unsigned long ant_pos_z
 )
{
  unsigned long X, Y, Z, ID;

  X = threadIdx.x + blockIdx.x * blockDim.x;
  Y = threadIdx.y + blockIdx.y * blockDim.y;
  Z = threadIdx.z + blockIdx.z * blockDim.z;
  ID = X + Y * Nx + Z * Nx * Ny;

  if(ant_pos_z-2 <= Z && Z <= ant_pos_z+2 &&
      X == ant_pos_x   && Y == ant_pos_y)
  {
    Ez[ID] = 0.0;
  }
  if(X == ant_pos_x && Y == ant_pos_y && Z == ant_pos_z)
  {
    Ez[ID] = wave;
  }
}
