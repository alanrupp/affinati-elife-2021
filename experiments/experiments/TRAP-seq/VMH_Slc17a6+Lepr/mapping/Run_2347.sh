## Aligning Vglut2/Lepr FASTQ to STAR genome
# 52 bp reads

FILES=`ls | grep Sample`

# fastq quality filter (q = 20)
for f in $FILES; do
  date
  echo 'Filtering' $f '('`echo $(($(echo ${FILES[@]/$f//} | cut -d/ -f1 | wc -w | tr -d ' ')+1))` 'of' `echo $FILES | wc -w`')'
  zcat ${f}/*.fastq.gz | fastq_quality_filter -q 20 -z -o ${f}/${f}_filtered.fastq.gz
done

# map reads to STAR
mkdir STAR_outs
cd STAR_outs

for f in $FILES; do
  date
  echo 'Mapping' $f '('`echo $(($(echo ${FILES[@]/$f//} | cut -d/ -f1 | wc -w | tr -d ' ')+1))` 'of' `echo $FILES | wc -w`')'
  ~/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir ~/STAR/GRCm38-92_Flpo_Cre_GfpL10a --sjdbGTFfile ~/STAR/GRCm38.92_Flpo_Cre_GfpL10a.gtf --readFilesCommand zcat --readFilesIn ../${f}/*_filtered.fastq.gz --outFileNamePrefix $f --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic
done
