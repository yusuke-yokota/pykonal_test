# pykonal_test
This code was modified by PyKonal, the Fast Marching Method (FMM; Sethian *et al.*, 1996) for solving the eikonal equation in Cartesian or spherical coordinates in 2 or 3 dimensions. The method implements mixed first- and second-order finite differences.
## Citation
If you make use of this code in published work, please cite White *et al.* (2020).

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
## Bugs
Please report bugs, feature requests, and questions through the [Issues](https://github.com/malcolmw/pykonal/issues "PyKonal Issues tracker") tracker.

## Py for GNSS-A
1. Cartesian_3Dfig3.py: 直交座標系の音響伝播計算
2. Spherical_3Dfig5.py: 球座標系の音響伝播計算（半径6371kmの球の表面5kmの球殻のみ計算される）
3. Path_fig1.py: 音響パスの計算（震源と受信点を拘束するのは難しくまだできていない）
4. VV.py: 傾斜パラメータが与えられるときの場を計算する（横田，2019，海洋音響学会誌）
5. mksvp.py: SVモデル作成用

## References
1. Sethian, J. A. (1996). A fast marching level set method for monotonically advancing fronts. *Proceedings of the National Academy of Sciences, 93*(4), 1591–1595. https://doi.org/10.1073/pnas.93.4.1591
2. White, M. C. A., Fang, H., Nakata, N., & Ben-Zion, Y. (2020). PyKonal: A Python Package for Solving the Eikonal Equation in Spherical and Cartesian Coordinates Using the Fast Marching Method. *Seismological Research Letters, 91*(4), 2378-2389. https://doi.org/10.1785/0220190318

