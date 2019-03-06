# WHAT?

三次元電磁界FDTD法による電場の数値解析プログラム
数値条件として、空間の中央に配置された1セルの高さを持つモノポールアンテナからガウスパルスを放射する。


## Compile

FDTD22法は-DSCHEME=2, FDTD24法は-DSCHEME=4をつけてコンパイル
```
cc main.c -fopenmp -std=c99 -O3 -DDEBUG -DSCHEME=2 -o FDTD22.exe
cc main.c -fopenmp -std=c99 -O3 -DDEBUG -DSCHEME=4 -o FDTD24.exe
```

## Usage

### Calculate

Nx,Ny,Nzはそれぞれ解析空間の幅、奥行き、高さをセルの個数で表す。
Ntは電磁界の更新回数を表す。
```
FDTD22.exe Nx Ny Nz Nt
```

1001ステップの計算を101x101x102の解析サイズで行う場合
```
FDTD22.exe 101 101 102 1001
```

### Visualize

出力される結果ファイルは`20190228042128_EM_FDTD22.bin`のように日時の文字列が先頭についたファイル名となっている。
中身はバイナリファイルであり、z = Nz/2となるxy断面におけるZ方向の電界の分布を単精度の浮動小数点で記録している。
`binToFigure.m`で実装されているbinToFigure関数により表面プロットを行う

```
binToFigure('20190228042128_EM_FDTD22.bin', 101, 101);
```
![1001ステップでの解析結果](https://github.com/NEETFUTURE/EM_FDTD/1001.png)

