# pykonal_test
This code was modified by PyKonal, the Fast Marching Method (FMM; Sethian, 1996 & White et al., 2020)  for calculate GNSS-A sound speed structure in order to solve the eikonal equation in Cartesian or spherical coordinates in 2 or 3 dimensions. The method implements mixed first- and second-order finite differences.
Using this code, we can calculate GNSS-A scale sound speed structure and GNSS-A data for simulation using GARPOS (Watanabe et al., 2020, https://zenodo.org/record/3992688#.X2xN4Xn7QuV).

## Citation
If you make use of this code in published work, please cite Sethian (1996) & White *et al.* (2020).
Citation of this simulation code itself is in preparation.

## Installation
https://malcolmw.github.io/pykonal-docs/installation.html

### Conda—recommended
```bash
sh$> conda create --name pykonal -c conda-forge python=3.8 'cython>=0.29.14' h5py numpy scipy
sh$> conda activate pykonal
(pykonal)sh$> cd path/to/pykonal
(pykonal)sh$> pip install .
```
### PIP
```bash
sh$> pip install cython>=0.29.14 h5py numpy scipy
sh$> cd path/to/pykonal
sh$> pip install .
```

## Py for GNSS-A
1. Cartesian_3Dfig3.py: 直交座標系の音響伝播計算
2. Spherical_3Dfig5.py: 球座標系の音響伝播計算（半径6371 kmの球の表面N kmの球殻のみ計算される）
3. Path_fig1.py: 音響パスの計算（震源と受信点を拘束するのは難しく,まだできていない）
4. VV.py: 傾斜パラメータが与えられるときの場(屈折非考慮)を計算する（横田，2019，海洋音響学会誌）
5. mksvp.py: SVモデル作成用
6. Simulator_*: GARPOSでのシミュレータ用，実行方法は下記

## Simulator
```bash
sh$> Simulator_3Dfig9.py
sh$> Simulator_addNoise.py
sh$> cd **/simulator
sh$> sh job.sh
sh$> (sh bat_mkfig.sh)
```

### Simulator_3D
1. 用意するファイル：svpC.csv 基本となる音速度ファイル（1500m/s均質で良い），Settings.ini セッティングファイル（詳しくは本体を見て下さい），linsample.csv 測線とデータ頻度を与える
2. 出てくるファイル：IMAG.1911.kaiyo_k4-obs.csv GARPOS用データファイル，IMAG.1911.kaiyo_k4-initcfg.ini GARPOS用cfgファイル

### Simulator_addNoise
1. 3D実施後同一階層で実行
2. 出てくるファイル：IMAG.${DAY}.kaiyo_k4-obs.csv ノイズ付与GARPOS用データファイル（ノイズ設定はSettings.ini），IMAG.${DAY}.kaiyo_k4-initcfg.ini GARPOS用cfgファイル

### job.sh & bat_mkfig.sh
1. GARPOS用のフォルダならどこでも良い．データとcfgを正しく設置する
2. demo.py (singleとfix)を通常通り回す > 初期設定ではdemo_prep/とdemo_res/に結果
3. bat_mkfig.shで作図

## References
1. Sethian, J. A. (1996). A fast marching level set method for monotonically advancing fronts. *Proceedings of the National Academy of Sciences, 93*(4), 1591–1595. https://doi.org/10.1073/pnas.93.4.1591
2. White, M. C. A., Fang, H., Nakata, N., & Ben-Zion, Y. (2020). PyKonal: A Python Package for Solving the Eikonal Equation in Spherical and Cartesian Coordinates Using the Fast Marching Method. *Seismological Research Letters, 91*(4), 2378-2389. https://doi.org/10.1785/0220190318

