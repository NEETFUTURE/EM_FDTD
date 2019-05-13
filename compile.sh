nvcc ./main.cu -o FDTD22.cuda -O2 -DSCHEME=2  -arch=sm_60 -lcuda,cudart -Xptxas -dlcm=ca -I /usr/local/cuda-8.0/samples/common/inc
nvcc ./main.cu -o FDTD24.cuda -O2 -DSCHEME=4  -arch=sm_60 -lcuda,cudart -Xptxas -dlcm=ca -I /usr/local/cuda-8.0/samples/common/inc
