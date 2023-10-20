#!/bin/bash 

for sample in Ribo_WT-1 Ribo_WT-2 Ribo_dMazEF-1 Ribo_dMazEF-2 

do 

	echo $sample 

	bamCoverage --bam "../mapped_bam/clean/"$sample"_clean.bam" -o $sample"_forward_coverage.bedgraph" --outFileFormat bedgraph --binSize 1 --filterRNAstrand forward
	
	bamCoverage --bam "../mapped_bam/clean/"$sample"_clean.bam" -o $sample"_reverse_coverage.bedgraph" --outFileFormat bedgraph --binSize 1 --filterRNAstrand reverse
	
	bamCoverage --bam "../mapped_bam/clean/"$sample"_clean.bam" -o $sample"_forward_coverage_Asite.bedgraph" --outFileFormat bedgraph --binSize 1 --Offset -12 --filterRNAstrand forward
	
	bamCoverage --bam "../mapped_bam/clean/"$sample"_clean.bam" -o $sample"_reverse_coverage_Asite.bedgraph" --outFileFormat bedgraph --binSize 1 --Offset -12 --filterRNAstrand reverse

done 

echo "done" 




