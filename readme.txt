
usage:
	GSVMining -option variable ...

Options:
	Required:
	-i/--analysisIDList	analysis id list in text file (indentical to id of raw read file)
	-c/--chromosomeList	chromosome list in txt file (per chromosome per line)
	-r/--readsDirectory	directory of raw read files (ex: id.read_1.fastq id.read_2.fastq)
	-f/--referencePath	path of reference in fasta format
	-o/--outputDirectory	directory of output
	Optional:
	-g/--gapThreshold	gap allowed between two blast alignment [default: 12]
	-d/--distThreshold	distance allowed between two blast alignment [default: 24]
	-e/--rpkmThreshold	rpkm threshold to filter low coverage candidates [default: 1]
	-s/--spanThreshold	span threshold to filter false positive candidates [default: 200]
	-l/--maxFragLength	maximal fragment length [default: 1000]
 
Note:	Pipeline for Genome Structural Variant Analysis
	1. samtools, bedtools, bwa, blast, perl should be already installed.
	2. blastdb should be built first.
 
	Fei Sang
	Deep Seq, Queen's Medical School
	School of Life Sciences, University of Nottingham
	Nottingham, Nottinghamshire, NG7 2JE, United Kingdom
	Email: fei.sang@nottingham.ac.uk
