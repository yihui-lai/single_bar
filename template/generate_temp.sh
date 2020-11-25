#!/bin/bash
mergefile="merge.sh"
main_direc=$PWD/..
zA="10"
for zlength in 100
do
hflth=`echo "${zlength}/2" | bc`
for zangle in 90 
#$(seq 75 10 165)
do
echo "length $zlength"
echo "half $hflth"
echo "angle $zangle"
special="BGO_${zangle}degree_${zlength}mm_A${zA}mm_sipm6mm"
Particle="muon"
Pa="mu-"
#Particle="elec"
#Pa="e-"
number="50"
for energy in 1
do
store_place="${main_direc}/${Particle}_${energy}GeV_N${number}_${special}/"
echo "make ${store_place}"
if [ ! -d "${store_place}" ]; then
   mkdir ${store_place}
fi
temp_name="${store_place}/template_${zangle}degree_${zlength}mm.cfg"
source temp_template.sh $zangle $zlength $zA 6 > ${temp_name}
ex_name="${store_place}/ex_${energy}GeV_N${number}_${Particle}_${zangle}degree_${zlength}mm.sh"
echo "hadd -f ${Particle}_${energy}GeV_N${number}_${special}.root ${Particle}_${energy}GeV_N${number}_${special}/*root" >> $mergefile
if [ -f ${ex_name} ]; then
   rm ${ex_name}
fi
for tag in $(seq 0 1 0)  
do
  mac_name="${store_place}run_${energy}GeV_N${number}_${Particle}_${tag}.mac"
  source temp_mac.sh $Pa $energy $number 1 -10> ${mac_name}
  source temp_jdl.sh ${main_direc} ${store_place} ${Particle} $energy $number $tag $temp_name ${mac_name} > ${store_place}condor_jobs_${energy}GeV_N${number}_${Particle}_${tag}.jdl
  echo "condor_submit ${store_place}condor_jobs_${energy}GeV_N${number}_${Particle}_${tag}.jdl" >> $ex_name
done
#source ${ex_name}
done
done
done
cd ..

cd ./template
