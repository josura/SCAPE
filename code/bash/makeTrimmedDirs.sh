shopt -s nullglob
for direc in */; do
   trimmed="trimmed"
   rm  $direc/*.trimmed.fastq
   echo mkdir "$direc$trimmed"
   mkdir "$direc$trimmed"
done

