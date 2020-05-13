# verify_bed

## What does this app do?
Code created to verify output bed files from nirvana_gff_to_bed app.</br></br>
Takes .gtf file of transcripts downloaded from UCSC table browser, reformats the file structure and then compares transcript coordinates against those in generated bed.

## What are typical use cases for this app?
Used for verifying bed files.

## What data are required for this app to run?
This app requires the generated bed file, and a .gtf file from UCSC table browser of uniq transcripts taken from the bed file.

## What does this app output?
It outputs 3 files:</br>
- ucsc_CDS_base_adjusted.gtf = input .gtf with required modifications for verifying, useful to sanity check changes
- only_in_ucsc.bed = bed file of transcripts only in UCSC .gtf file (i.e. non matching)
- ucsc_nirvana_match_transcripts.bed = bed file of matching transcripts 
