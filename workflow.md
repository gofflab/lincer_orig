# Offline
1) Assemble transcriptomes (.gtfs) for each of your samples
2) Create merged assembly (cuffmerge or stringtie-merge)
3) Create 'mask' .gtf file to exclude ANY gene with existing annotation that is NOT a lincRNA (starting from the most recent reference .gtf file)

# From the top
1) Find a way to estimate max coverage for each transcript (not gene, not exon) across all your samples
2) Same for FPKM if possible (can use as a filter later)
3) Filter for:
	- min 3x coverage
	- >250bp
	- >1 exon

4) take "putative lncRNAs" going forward

A) If you care whether or not they are lncRNAs, do the following:


	5) PFam (retain inforamtion for later use)
	6) PhyloCSF
		- hard to run directly, look into UCSC trackHub (ask Loyal to dig for it if you can't find it) to scrape largest lod score
		-
	7) Flag as 'pseudogene' if overlapping with annotated pseudogene

B) otherwise, do your differential analysis first!  Then only go forward with those genes that have meaningfull differential expression.
	- Then take them forward through lincRNA pipeline. 

