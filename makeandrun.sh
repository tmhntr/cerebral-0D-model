#!/bin/bash
make veryclean
make cbf

D=$(date "+%F")

if [[ ! -d ../../outputs/$D ]]
then
  mkdir ../../outputs/$D
fi

n=0
fol=$(printf "run%02d\n" $n)
while [[ -d ../../outputs/${D}/${fol} ]]
do
  n=$(($n + 1))
  fol=$(printf "run%02d" $n)
done

mkdir ../../outputs/${D}/${fol}
mv cbf ../../outputs/${D}/${fol}
cd ../../outputs/${D}/${fol}

./cbf $1 $2 $3 $4 $5
