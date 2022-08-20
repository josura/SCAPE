shopt -s nullglob
for direc in */; do
   aligned="alignedToTranscript"
   #rm  $direc/*.sam
   echo mkdir "$direc$aligned"
   mkdir "$direc$aligned"
done

