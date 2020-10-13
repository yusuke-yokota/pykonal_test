#!/bin/bash -eu
#./Simulator_3Dfig9.py
./Simulator_3D_v2.py
./Simulator_addNoise2b.py
./Simulator_addNoise2c.py
cd ~/sgobs/garpos/simulator
./job.sh
cd -
