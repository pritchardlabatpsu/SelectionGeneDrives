
### Analyze sgRNAs using PoolQ

## Setup
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Files
idx.reads = "PC9OsimertinibScreenIndexReads.csv"

### Place sequencing files in a folder named SequencingFiles ###

## Prep sequencing files for PoolQ

# Unzip files - may need to run twice if OneDrive is being weird
system("gunzip SequencingFiles/*.fastq.gz")

# Read in index read file
idx.df = read.csv(idx.reads,header=T,stringsAsFactors=F)

# Generate pseudo barcode files - this takes a few minutes

seq.files = list.files(path="./SequencingFiles/",pattern="\\.fastq$")
system("mkdir -p ./IndexFiles/")

# May require that homebrew be added to PATH
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/opt/homebrew/bin/", sep=":")) # required to add gsed command to path; in future, find path with terminal command "which gsed"

for (i in 1:nrow(idx.df)) {
  
  sample.i = idx.df$Sample[i]
  brcd.i = idx.df$Index.Read[i]
  
  seq.file.i = grep(sample.i,seq.files,value=T) # identify appropriate sequencing file
  
  # Confirm that barcode in sequencing file matches index assignment
  
  
  # Generate pseudo index files
  system(paste0("gsed '2~4s/.*/",brcd.i,"/' ./SequencingFiles/",seq.file.i," > ./IndexFiles/",sample.i,"_idxtmp.fastq")) # First, replace sequencing reads with appropriate barcode
  system(paste0("gsed '4~4s/.*/JJJJJJJJ/' ./IndexFiles/",sample.i,"_idxtmp.fastq > ./IndexFiles/",sample.i,"_idx.fastq"))  # Then, replace Phred scores with string of J's
  system(paste0("rm ./IndexFiles/",sample.i,"_idxtmp.fastq")) # Remove temporary index file
}

# Concatenate sequence files and index files

system("mkdir -p ./PoolQResults/")
system(paste0("cat ./SequencingFiles/*.fastq > ./PoolQResults/concatenated_sequences.fastq"))
system(paste0("cat ./IndexFiles/*.fastq > ./PoolQResults/concatenated_idx.fastq"))

## Prep PoolQ input files

# Generate row reference file

# Import Brunello file
guides = read.csv("BrunelloGuides.tsv",sep="\t",header=T,stringsAsFactors=F)

# Reformat
guides.frmt = data.frame(Barcode = guides$sgRNA.Target.Sequence,
                         ID = guides$Target.Gene.Symbol)

# Save file
write.table(guides.frmt,"./PoolQResults/BrunelloReferenceFile.csv",row.names=F,col.names=F,sep=",")

# Generate column reference file

# Import index barcode file

idx.rfrmt = idx.df[,c("Index.Read","Sample")]
colnames(idx.rfrmt) = c("Barcodes","Condition")

write.table(idx.rfrmt,"./PoolQResults/ColReferenceFile.txt",row.names=F,col.names=F,sep="\t")

## Run PoolQ

setwd("./PoolQResults/")

system("../PoolQ/poolq3.sh --row-reads concatenated_sequences.fastq --col-reads concatenated_idx.fastq --col-reference ColReferenceFile.txt --row-reference BrunelloReferenceFile.csv --row-barcode-policy PREFIX:CACCG@16-31 --col-barcode-policy FIXED:0")

## Next, calculated LFC and p.values using MAGeCKFlute in MAGeCKFluteScript_PC9Osimertinib1.R
