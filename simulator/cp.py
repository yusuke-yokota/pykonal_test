#!/home/sgo/anaconda3/bin/python3
import os
import sys
import glob
import configparser
from optparse import OptionParser
import os
import sys

basfile='~/sgobs/garpos/simulator/initcfg/IMAG/IMAG.1911.kaiyo_k4-initcfg.ini'
for epi in range(100):
	ix = 2001 + epi % 12 + 100 * (epi // 12)
	cfgfile='~/sgobs/garpos/simulator/initcfg/IMAG/IMAG.2%s.kaiyo_k4-initcfg.ini' % ix
	os.system('sed -i -e "1i ' + hd + '" ' + of )
	
