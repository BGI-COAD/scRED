inf='bam_2.list'
outdir='/hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/mice_reprogram/2.RNA_edit_site/1.run_pipe_run2'
rm -f $outdir/qsub.sh
for line in `cat $inf`
do
prefix=`echo $line|awk -F "," '{print $1}'`
bam=`echo $line|awk -F "," '{print $2}'`
mkdir -p $outdir/$prefix
echo "samtools index $bam
sh /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/mice_reprogram/pipe_test/full_pipe_run/pipe.sh $bam /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/mice_reprogram/test/hard_ft/bin $outdir/$prefix /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/ref/mm10_tophat2/mm10.fa /hwfssz1/ST_PRECISION/USER/lvtianhang/software/bwa-0.6.2/bwa /zfssz3/ST_CANCER/CGR/SHARE/RNA_editing_pipeline/hard_filter/database/annovar /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/pipe/ref/simpleRepeat.bed /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/pipe/ref/mm10.alu.bed /hwfssz1/ST_PRECISION/USER/lvtianhang/software/samtools-1.9/bin/bin /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/mice_reprogram/ref/mm10_mrna.fa /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/mice_reprogram/ref/run/dbsnp.txt /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/mice_reprogram/ref/mm10_mrna.fa.fai /hwfssz1/ST_PRECISION/USER/lvtianhang/work/RNA_editing/mice_reprogram/ref/refGene.txt">$outdir/$prefix/${prefix}.sh
echo "qsub -wd ${outdir}/${prefix} -l vf=2g,num_proc=1 -P P18Z10200N0350 -binding linear:1 -q st.q ${outdir}/$prefix/${prefix}.sh">>${outdir}/qsub.sh  
done
