#!/bin/bash -eu
# -e エラーがあったらシェルスクリプトをそこで打ち止めにしてくれる
# -u 未定義の変数を使おうとしたときに打ち止めにしてくれる

sitelist="IMAG"

singledir="single-bsp331"
cfgfixdir="cfgfix-bsp331"
resultdir="result-bsp333"

finiprep="Settings-prep.ini"
finifix="Settings-fix.ini"
fininoise="Settings-noise.ini"

isingle=true
irun=true
irunfix=true

for site in ${sitelist}
do
	f1911="initcfg/${site}/${site}.1911.kaiyo_k4-initcfg.ini"
	for obs in `ls ./obsdata/${site}/${site}.[2-9]0*.*yo*obs.csv`
#	for f1 in `ls ./initcfg/${site}/${site}.[2-9][0-9]*.*yo*initcfg.ini`
	do
		bn=`basename $obs`
		f1="initcfg/${site}/${bn/obs.csv}initcfg.ini"
#		obs="obsdata/${site}/${bn/initcfg.ini}obs.csv"
		./demo.py -i $finiprep -f $f1 -d demo_prep/$site --maxcore 1
	done
#	cp demo_prep/IMAG/IMAG.1911.kaiyo_k4-res.dat demo_prep/IMAG/IMAG.1910.kaiyo_k4-res.dat
	./makeFixCfg.py -d cfgfix --res_singles "demo_prep/${site}/*res.dat"
	for f2 in `ls ./initcfg/${site}/${site}.[2-9]0*.*yo*initcfg.ini`
	do
		bn=`basename $f2`
		cfgfix_ini="cfgfix/${site}/${bn/initcfg.ini}fix.ini"
		./demo.py -i $finifix -f $cfgfix_ini -d demo_res/$site --maxcore 2
	done
	
done
./bat_mkfig.sh
