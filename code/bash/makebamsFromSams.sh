shopt -s nullglob
for direc in */; do
   aligned="aligned/result.sam"
   echo samtools view -S -b "$direc$aligned" > "${direc}aligned/result.bam"
   samtools view -S -b "$direc$aligned" > "${direc}aligned/result.bam"
done

