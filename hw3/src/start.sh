rm *.out *.err
make
sbatch job-upc
sleep 15
clear
tail -10 upc-HW3.*.out
