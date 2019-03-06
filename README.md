# EM_FDTD

三次元電磁界FDTD法プログラム
数値条件として、空間の中央から1セルの高さを持つモノポールアンテナから、ガウスパルスを放射する。


## Compile

FDTD22法は-DSCHEME=2, FDTD24法は-DSCHEME=4をつけてコンパイル
```
cc main.c -fopenmp -std=c99 -O3 -DDEBUG -DSCHEME=2 -o FDTD22.exe
cc main.c -fopenmp -std=c99 -O3 -DDEBUG -DSCHEME=4 -o FDTD24.exe
```

## USAGE

Nx,Ny,Nzはそれぞれ解析空間の幅、奥行き、高さをセルの個数で表す。
Ntは電磁界の更新回数

```
FDTD22.exe Nx Ny Nz Nt
```

1001ステップの計算を101x101x102の解析サイズで行う場合
```
FDTD22.exe 101 101 102 1001
```
