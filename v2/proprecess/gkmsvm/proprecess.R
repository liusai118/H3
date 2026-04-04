###############
###############
###############
library(GenomicFeatures)
library(Homo.sapiens)
library(clusterProfiler)
library(dplyr)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9.masked)
library(gkmSVM)
library(data.table)
pick_by_label_with_gap <- function(df, gap = 100L) {
  stopifnot(all(c("seqnames","start","end","label") %in% names(df)))
  gap <- as.integer(gap)
  
  df <- df %>%
    mutate(
      seqnames = as.character(seqnames),
      start = as.integer(start),
      end   = as.integer(end),
      label = as.numeric(label)
    ) %>%
    filter(!is.na(seqnames), !is.na(start), !is.na(end), !is.na(label)) %>%
    filter(start <= end)
  
  
  res <- df %>%
    group_by(seqnames) %>%
    group_modify(~{
      x <- .x %>% arrange(desc(label), start, end)
      

      ir <- IRanges(start = x$start, end = x$end)
      

      pad <- gap - 1L
      ir_pad <- IRanges(start = pmax(1L, start(ir) - pad),
                        end   = end(ir) + pad)
      
      keep <- logical(nrow(x))
      selected_pad <- IRanges()  
      
      for (i in seq_len(nrow(x))) {
        if (length(selected_pad) == 0L) {
          keep[i] <- TRUE
          selected_pad <- c(selected_pad, ir_pad[i])
        } else {

          if (countOverlaps(ir_pad[i], selected_pad) == 0L) {
            keep[i] <- TRUE
            selected_pad <- c(selected_pad, ir_pad[i])
          }
        }
      }
      
      x[keep, , drop = FALSE]
    }) %>%
    ungroup()
  
  res
}


check_gap_ok <- function(df, gap = 100L) {
  df <- df %>% arrange(seqnames, start, end)
  df %>%
    group_by(seqnames) %>%
    summarize(
      min_gap = {
        if (n() < 2) NA_integer_
        else min(start[-1] - end[-n()])
      },
      ok = {
        if (n() < 2) TRUE
        else all(start[-1] - end[-n()] >= gap)
      },
      .groups = "drop"
    )
}

write_sequences_to_fasta <- function(sequences, bed_file, filepath) {
  stopifnot(length(sequences) == nrow(bed_file))
  
  con <- file(filepath, open = "w")
  on.exit(close(con), add = TRUE)
  
  headers <- paste0(">", bed_file$seqnames, ":", bed_file$start, "-", bed_file$end)
  seqs <- unname(as.character(sequences))
  seqs <- gsub("\\s+", "", seqs)
  

  fasta_lines <- as.vector(rbind(headers, seqs))
  writeLines(fasta_lines, con)
}

###############
###############
###############
###############
genome <- BSgenome.Mmusculus.UCSC.mm9
a <- fread("H3_minOv2_summits250_consensusPeaks.bed")
colnames(a)
a[, mid := (V2 + V3) %/% 2]
a$start <- a$mid -49
a$end <- a$mid +50
a <- a[,c(1,8,9)]
colnames(a) <-  c("seqnames","start","end")
write.table(a,"H1.bed",col.names = F,row.names = F,sep = "\t",quote = F)
genNullSeqs('H1.bed',nMaxTrials=40,xfold=3,genome=BSgenome.Mmusculus.UCSC.mm9.masked,outputPosFastaFN='n1.fa', outputBedFN='neg1x_n3.bed', outputNegFastaFN='neg1x_n3.fa')
###############
###############
databed1 <- read.table("H1.bed")
bed_file2 <- import("neg1x_n3.bed", format = "BED")
databed2 <- read.table("neg1x_n3.bed")
colnames(databed2) <- c("seqnames","start","end")
genome <- BSgenome.Mmusculus.UCSC.mm9
sequences2 <- getSeq(genome, bed_file2)
###############
pos <- read.table("H1.bed")[,1:3]
pos <- pos[!duplicated(pos$V2),]
neg <- read.table("neg1x_n3.bed")
neg <- neg[!duplicated(neg$V2),]
pos$label <- 1
neg$label <- 0
all <- rbind(pos,neg)
all <- all[!duplicated(all$V2),]

colnames(all) <- c("seqnames","start","end","label")
############### 
library(dplyr)
all2 <- all %>%
  filter(!grepl("[eE]", as.character(end)))
all_unique <- pick_by_label_with_gap(all2, gap = 100L)


pos1 <- all_unique %>% filter(label == 1)
all_unique <- all_unique[order(all_unique$label,decreasing = T),]


pos1 <- all_unique[all_unique$label == 1,]
neg <-  all_unique[all_unique$label == 0,]
neg$start <- neg$start -1 
neg <- neg[!duplicated(neg$start),]
neg1 <- neg[sample((nrow(neg)), 1*(length(which(all_unique$label == 1)))), ]
class(all2$end)
all2 <- rbind(pos1,neg1)
nrow(all2)

all2 <- all2[order(all2$label,decreasing = T),]
all2$end <- all2$end - 1 
write.table(all2,"all.bed",col.names = F,row.names = F,sep = "\t",quote = F)
bed_file3 <- import("all.bed", format = "BED")
library(BSgenome.Mmusculus.UCSC.mm9)
genome <- BSgenome.Mmusculus.UCSC.mm9
sequences3 <- getSeq(genome, bed_file3)
all2_bed <- all2[,1:3]
write_sequences_to_fasta(sequences3,all2_bed, filepath = "all1.fasta")
all<- read.table("all1.fasta")
fastq <- readLines("all1.fasta")
seq_lines <- fastq[seq_along(fastq) %% 2 == 0]
all2$sequence <- seq_lines
write.csv(all2,"all.csv")

