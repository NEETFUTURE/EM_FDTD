for i in `seq 5`;do
./FDTD22.cuda 128 128 128 1000
./FDTD24.cuda 128 128 128 1000
done
