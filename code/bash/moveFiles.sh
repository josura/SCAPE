shopt -s nullglob
for f in *_2.fastq; do
   fshorted="${f%_2.*}"
   mv "$f" "$fshorted/$dest"
   echo mv "$f" "$fshorted/$dest"
done
