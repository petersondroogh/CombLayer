#! /bin/bash
#
# first argument is the xml file with geometry variables

if [ $# != 1 -a $# != 2 ]; then
        echo "Usage: run.sh geometry.xml [TW]"
        exit 2
fi

xml=$1
if [ ! -e $xml ]; then
    echo "xml file $xml does not exist"
    exit 1
fi
pstudy-nowrap -i $xml -setup -job flat.xml

TW=""
if [ $2 = "TW" ]; then
    TW="-TW -WP 0 0 -11.45"
fi

for dir in case*
do
    cd $dir
     x5=$(getvariable flat.xml  F5X);   y5=$(getvariable flat.xml  F5Y);  z5=$(getvariable flat.xml  F5Z)
    x15=$(getvariable flat.xml F15X);  y15=$(getvariable flat.xml F15Y); z15=$(getvariable flat.xml F15Z)
    x25=$(getvariable flat.xml F25X);  y25=$(getvariable flat.xml F25Y); z25=$(getvariable flat.xml F25Z)

    x35=$(getvariable flat.xml F35X);  y35=$(getvariable flat.xml F35Y); z35=$(getvariable flat.xml F35Z)
    x45=$(getvariable flat.xml F45X);  y45=$(getvariable flat.xml F45Y); z45=$(getvariable flat.xml F45Z)
    x55=$(getvariable flat.xml F55X);  y55=$(getvariable flat.xml F55Y); z55=$(getvariable flat.xml F55Z)


    test $(hostname) \!= "jasper" && module load gcc/4.9.2


    ess -r -nF5 6 -lowMod Butterfly -topMod Butterfly -lowModFlowGuide On -topModFlowGuide On -lowWaterDisc On -topWaterDisc On  -topPreCooling None  -x flat.xml -T point free $x5 $y5 $z5 -T point free $x15 $y15 $z15 -T point free $x25 $y25 $z25 -T point free $x35 $y35 $z35 -T point free $x45 $y45 $z45 -T point free $x55 $y55 $z55  -T surface object Bulk 3 -TMod energy 0 "5E-9 20E-9 100E-9" -v TopPreMaterial2 4005 --validCheck 1000 $TW butterfly


    test $(hostname) \!= "jasper" && module unload gcc/4.9.2
    sed -i -n '/^sdef/{p;:a;N;/prdmp/!ba;s/.*\n/si1 h -7 7\nsp1 0 1\nsi2 h -1.6 1.6\nsp2 0 1\n/};p' butterfly1.x
    remove-weights butterfly1.x

    mv butterfly1.x inp
    sed -i -e "s/nps 20000/stop f5 0.003/" inp
#    sed -i -e "s/nps 20000/stop nps 10/" inp
    sed -i -e "s/^t/c t/" inp

    setomega inp F5 5
    setomega inp F15 15
    setomega inp F25 25

    setomega inp F35 35
    setomega inp F45 45
    setomega inp F55 55
#    rm -f ObjectRegister.txt  Renumber.txt  Spectrum.log
    remove-weights inp
    cd ..;
done
