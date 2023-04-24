

outdir=/rds/general/user/ah3918/ephemeral/PD_PILOT_MULTIOME/roche_data

for i in {1..4}
do
  cp $outdir/PD00${i}/outs/filtered_feature_bc_matrix.h5 $outdir/Colossus/PD${i}_filtered_feature_bc_matrix.h5
  cp $outdir/PD00${i}/outs/atac_fragments.tsv.gz.tbi $outdir/Colossus/PD${i}_atac_fragments.tsv.gz.tbi
  cp $outdir/PD00${i}/outs/atac_fragments.tsv.gz $outdir/Colossus/PD${i}_atac_fragments.tsv.gz
done

 ‘GenomeInfoDb’, ‘GenomicRanges’, ‘IRanges’, ‘Rsamtools’, ‘S4Vectors’, ‘BiocGenerics’
