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
  fol=$(printf "prodrun%02d" $n)
done

mkdir ../../outputs/${D}/${fol}
mv cbf ../../outputs/${D}/${fol}
cd ../../outputs/${D}/${fol}

run=1
for cow in 0 1 2 3 4 5
do
  for HR in 50 70 90
  do
    for AF in 0 1
    do
      echo "./cbf $run $AF $cow $HR 1" >> job_file
      run=$(($run + 1))
    done
  done
done

FILE="cerebral_run${n}.sh"
# the EOM has to be first character on line, else file is assumed not to be finished. this comment does go into the file, so moving it elsewhere.
/bin/cat <<EOM >$FILE
#!/bin/bash
#SBATCH --account=def-kharches
#SBATCH --mail-user=thunte27@uwo.ca		# Email
#SBATCH --mail-type=ALL								# Email me notifications
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32 			# Graham has 32 cores, Cedar has 48 cores to a node.
#SBATCH --mem-per-cpu=1024M      		# memory; default unit is megabytes
#SBATCH --time=00-00:59          			# time (DD-HH:MM)
parallel -j 32 < job_file
EOM
sbatch $FILE
