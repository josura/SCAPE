shopt -s nullglob
for direc in */; do
   aligned="aligned"
   #rm  $direc/*.sam
   echo ls "$direc$aligned"
   ls "$direc$aligned"
done

