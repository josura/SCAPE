echo "sample, quant_file, condition" >> sample.csv
shopt -s nullglob
for direc in */; do
   quantfile="${direc}salmon/quantification/quant.sf"
   sample=${direc%/}
   #rm  $direc/*.sam
   echo "${sample},$(pwd)/${quantfile},BoneMarrowALL" >> sample.csv
done
