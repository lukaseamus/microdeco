#clean feild trial 16S
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

path <- "/Users/gengen/Documents/nurdi/nurdi_seqs"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

FWD <- "CCTAYGGGRBGCASCAG" 
REV <- "GGACTACNNGGGTATCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, matchIDs=TRUE, maxN = 0, multithread = TRUE, verbose = T)
#This step took 1.75hrs for 64GB of data

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
#check hits for 10th sample in list
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

t1=Sys.time()
cutadapt <-"/opt/anaconda3/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
t2=Sys.time()
t2-t1
#3.6 hrs for 90GB

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
plotQualityProfile(cutFs[3])
plotQualityProfile(cutRs[3])

t1=Sys.time()
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out=filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                  truncLen=c(260,212), 
                  minLen=200, 
                  maxN=0, 
                  truncQ=2, 
                  maxEE=c(4,4),
                  rm.phix=TRUE, 
                  compress=TRUE, 
                  verbose=TRUE, 
                  multithread=T)
head(out)
t2=Sys.time()
t2-t1
write.csv(out,file="/Users/gengen/Documents/nurdi/nurdi_filterntrim.csv")
#lost 50% of reads so changing maxEE to (4,4). This retained on average 77% of reads (range = 60%-95%)

t1=Sys.time()
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
t2=Sys.time()
t2-t1

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


t1=Sys.time()
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
t2=Sys.time()
t2-t1



#up to here
t1=Sys.time()
dadaFs <- dada(derepFs, err = errF, multithread = TRUE,verbose = T)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE,verbose = T)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, "/Users/gengen/Documents/nurdi/nurdi_seqtab.rds") 
t2=Sys.time()
t2-t1

t3=Sys.time()
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
t4=Sys.time()
t4-t3
table(nchar(getSequences(seqtab.nochim)))
saveRDS(seqtab.nochim, "/Users/gengen/Documents/nurdi/nurdi_seqtab.nochim.rds") 

taxa <- assignTaxonomy(seqtab.nochim, "/Users/gengen/Documents/Taxonomy_databases/SILVA/silva_nr_v132_train_set.fa.gz",multithread=TRUE)
#started at 4:36 AM - took 10 mins
head(taxa)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#save important files
write.csv(track,file="/Users/gengen/Documents/nurdi/nurdi_reads_thru_pipeline.csv")
saveRDS(taxa, "/Users/gengen/Documents/nurdi/nurdi_tax.rds") # CHANGE ME to where 





