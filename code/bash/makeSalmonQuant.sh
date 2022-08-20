shopt -s nullglob
for direc in */; do
   trimmed="salmon"
   toconvert="aligned/result.sam"
   result="salmon/result.fastq"
   echo mkdir "$direc$trimmed"
   #mkdir "$direc$trimmed"
   echo "samtools fastq $direc$toconvert > $direc$result"
   #samtools fastq $direc$toconvert > $direc$result
   echo "salmon quant -i ../../referenceGenomes/human/referenceTranscriptome/partialIndexSalmon/default/ -l A -1 ${direc}trimmed/reads_1.trimmed.fastq -2 ${direc}trimmed/reads_2.trimmed.fastq --validateMappings -o ${direc}salmon/quantification"
   #salmon quant -i ../../referenceGenomes/human/referenceTranscriptome/partialIndexSalmon/default/ -l A -1 ${direc}trimmed/reads_1.trimmed.fastq -2 ${direc}trimmed/reads_2.trimmed.fastq --validateMappings -o ${direc}salmon/quantification
   salmon quant -i ../../referenceGenomes/human/referenceTranscriptome/fullIndexSalmon/default/ -l A -1 ${direc}trimmed/reads_1.trimmed.fastq -2 ${direc}trimmed/reads_2.trimmed.fastq --validateMappings -o ${direc}salmon/quantificationFull
done

