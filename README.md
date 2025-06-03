# SunflowerMutagenesisAssembly

High level summary:
Genome Assembly - HiFi + OmniC:
*Software used in this pipeline was often installed via conda to an environment, then loaded with "source /mmfs1/home/USER.NAME/miniconda3/bin/activate envNameHere"
1. Received raw data - 2 PacBio Revio HiFi SMRT cells + R1/R2 OmniC
2. Merge/concatenate 2 SMRT cell data into one fq.gz file.
3. Use HifiAdaptFilt to ensure that no adapters/contaminants remain (they didn't).
	bash hifiadapterfilt.sh -p mergedHiFiReads.fq.gz -t #threads
4. Use meryl then genomescope2.0 to visualize k-mer counts and get homo/het, haplotype, repeat, and error stats
	meryl-1.4.1/bin/meryl count k=21_OR_31 threads=#threads memory=#gbRAM output hifi_kmers_k21_OR_31.meryl mergedHiFiReads.fq.gz
	Rscript genomescope2.0/genomescope.R -i input_merylHistogram_file -o output_dir -p ploidy/2 -k kmer_length/31_OR_21
		Then look at the resulting plots and summary.txt file!
		*My understanding is k21 is standard/common especially for Illumnina, but VGP pipeline suggested 31 for HiFi. We tried both.
5. Run hifiasm using merged hifi reads (fq.gz) and both R1/R2 omniC (fq.gz)
	hifiasm -o outputLabel.asm -t #threads --h1 hic/omnicR1.fq.gz --h2 hic/omnicR2.fq.gz mergedHiFiReads.fq.gz
		*p_ctg means polished contigs, and is output for primary, hap1 and hap2 when in HiC-phased mode
6. Run gfastats to get contig statistics and convert .gfa to .fasta format
	gfastats/build/bin/gfastats primary_OR_hap1_OR_hap2.p_ctg.gfa --threads #threads -o primary_OR_hap1_OR_hap2.p_ctg.fasta
		Then look at the statistics in the accompanying output - #contigs/scaffolds, length, N50, etc.
7. Run BUSCO on contigs
	busco -i primary_OR_hap1_OR_hap2.p_ctg.fasta -m genome -c #threads -o outputDirectory -l eudicots_odb10
		*eudicots is the closest lineage to sunflower for BUSCO gene completeness analysis
8. Run Merqury on contigs
	merqury/merqury.sh hifi_kmers_k21_OR_31.meryl primary_OR_hap1_OR_hap2.p_ctg.fasta outputDirectory
9. Generate contact map for contigs
	# Index the fasta contigs file
	samtools faidx primary_OR_hap1_OR_hap2.p_ctg.fasta
	# Print first two columns of index into genome file
	cut -f1,2 primary_OR_hap1_OR_hap2.p_ctg.fasta.fai > primary_OR_hap1_OR_hap2.p_ctg.genome
	# Index the fasta
	bwa-mem2/bwa-mem2 index primary_OR_hap1_OR_hap2.p_ctg.fasta

	# Align fasta to OmniC data
	bwa-mem2/bwa-mem2 mem -5SP -T0 -t#threads primary_OR_hap1_OR_hap2.p_ctg.fasta OmniCR1.fq.gz OmniCR2.fq.gz o primary_OR_hap1_OR_hap2_alignedToOmniC.sam
	# Find ligation junctions/parse
	pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 -nproc-in #threads -nproc-out #threads -chroms-path primary_OR_hap1_OR_hap2.p_ctg.genome primary_OR_hap1_OR_hap2_alignedToOmniC.sam > primary_OR_hap1_OR_hap2_alignedToOmniC_parsed.pairsam
	# Sort the parsed pairs
	pairtools sort --tmpdir=/path/to/tempdir --nproc #threads primary_OR_hap1_OR_hap2_alignedToOmniC_parsed.pairsam > primary_OR_hap1_OR_hap2_alignedToOmniC_sorted.pairsam
	# Mark duplicates
	pairtools dedup --nproc-in #threads --nproc-out #threads --mark-dups --output-stats primary_OR_hap1_OR_hap2_alignedToOmniC_pairtools_stats.txt --output primary_OR_hap1_OR_hap2_alignedToOmniC_dedup.pairsam primary_OR_hap1_OR_hap2_alignedToOmniC_sorted.pairsam
	# Split pairsam into two files
	pairtools split --nproc-in #threads --nproc-out #threads --output-pairs primary_OR_hap1_OR_hap2_alignedToOmniC_mapped.pairs --output-sam primary_OR_hap1_OR_hap2_alignedToOmniC_unsorted.bam primary_OR_hap1_OR_hap2_alignedToOmniC_dedup.pairsam 
	# Generate the final bam file
	samtools sort -@#threads -T /path/to/tempdir/temp primary_OR_hap1_OR_hap2_mappedToOmniC.PT.bam primary_OR_hap1_OR_hap2_alignedToOmniC_unsorted.bam 
	# Index the final bam file
	samtools index primary_OR_hap1_OR_hap2_mappedToOmniC.PT.bam

	# Generate contact map data
	samtools view -h primary_OR_hap1_OR_hap2_mappedToOmniC.PT.bam | PretextMap -o primary_OR_hap1_OR_hap2_omnic_pretext --sortby nosort --mapq 10
	# Visualize contact map
	PrextextSnapshot -m primary_OR_hap1_OR_hap2_omnic_pretext -f "png"
10. Scaffold with YaHS
	yahs/yahs primary_OR_hap1_OR_hap2.p_ctg.fasta primary_OR_hap1_OR_hap2_mappedToOmniC.PT.bam 
		*Purging duplicates isn't necessary before this because of the phasing with OmniC
11. Run BUSCO on scaffolds
	Just like above, but with scaffold fasta
12. Run Merqury on scaffolds
	Just like above, but with scaffold fasta
13. Generate contact map for scaffolds
	Just like above, but with scaffold fasta
14. Use JupiterPlot to compare scaffolds to a published reference genome
	The first 17 scaffolds should be MUCH larger than the rest, corresponding to chromosomes. Use them as input to JupiterPlot:
	jupiter name=OmniCprimary_OR_hap1_OR_hap2 ref=reference_17scaffs.fasta fa=primary_OR_hap1_OR_hap2.fasta t=#threads m=0 ng=0 maxScaff=17 labels=both
		*m is minimum reference chromosome size. Could use this instead of trimming ref fasta, e.g. 10-100Mb. ng limits to percent of ref genome size to map, disable this. maxScaff limits number of scaffs, can use 17 to exclude non-Chr ones. labels=both applies labels to both ref and query fasta
		*Do NOT use the resulting jpg. It was much lower quality than svg.
15. Relabel chromsomes according to links with reference. Chr01-17, Chr00c000 for scaffolds, moving Scaff18 to Chr00c001.
