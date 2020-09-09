#!/bin/bash -eu
# -e エラーがあったらシェルスクリプトをそこで打ち止めにしてくれる
# -u 未定義の変数を使おうとしたときに打ち止めにしてくれる

sitelist="IMAG"

resultdir="demo_res"
#resultdir="demo_prep"

#fini="Settings-single.ini"
#finifix="Settings-fix.ini"

for site in ${sitelist}
do
	for f in `ls ${resultdir}/${site}/${site}.1911*-res.dat`
#	for f in `ls ${resultdir}/${site}/${site}.[0-1][0-9]*-res.dat`
	do
		cmd="p06_plotResults_ts.py --site ${site} --resfile ${f}"
		echo ${cmd}
		${cmd}
		cmd="p07_Vd+iVs_dV.py --site ${site} --resfile ${f}"
		echo ${cmd}
		${cmd}
		if [ $? = 1 ]; then
			echo "error!"
			exit
		fi
	done
	sudo cp ${resultdir}/${site}/fig/*.*t.s*.png /mnt/owncloud_webdav/share/resiplot/${site}/fix
#	sudo cp ${resultdir}/${site}/fig/*.*t.s*.png /mnt/owncloud_webdav/share/resiplot/${site}/nofix
	sudo cp ${resultdir}/${site}/fig/*.*dvcor*.png /mnt/owncloud_webdav/share/resiplot/${site}/fix
#	sudo cp ${resultdir}/${site}/fig/*.*dvcor*.png /mnt/owncloud_webdav/share/resiplot/${site}/nofix
done

#resultdir="demo_res"
resultdir="demo_prep"

#fini="Settings-single.ini"
#finifix="Settings-fix.ini"

for site in ${sitelist}
do
	for f in `ls ${resultdir}/${site}/${site}.1911*-res.dat`
#	for f in `ls ${resultdir}/${site}/${site}.[0-1][0-9]*-res.dat`
	do
		cmd="p06_plotResults_ts.py --site ${site} --resfile ${f}"
		echo ${cmd}
		${cmd}
		cmd="p07_Vd+iVs_dV.py --site ${site} --resfile ${f}"
		echo ${cmd}
		${cmd}
		if [ $? = 1 ]; then
			echo "error!"
			exit
		fi
	done
#	sudo cp ${resultdir}/${site}/fig/*.*t.s*.png /mnt/owncloud_webdav/share/resiplot/${site}/fix
	sudo cp ${resultdir}/${site}/fig/*.*t.s*.png /mnt/owncloud_webdav/share/resiplot/${site}/nofix
#	sudo cp ${resultdir}/${site}/fig/*.*dvcor*.png /mnt/owncloud_webdav/share/resiplot/${site}/fix
	sudo cp ${resultdir}/${site}/fig/*.*dvcor*.png /mnt/owncloud_webdav/share/resiplot/${site}/nofix
done
