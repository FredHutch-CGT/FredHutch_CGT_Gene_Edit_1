#!/usr/bin/env Rscript
# Filename: kevin_genediting_mutations.R
#
# WARNINGS: This script was adapted from J Davison's original
#   Be wary of two coding styles intersecting. Things get weird.
#
# Targeted mutation analysis for HKiem lab.
# ==> Supports stitched, filtered, needle-aligned reads.
#
# J Davison 18mar2015/02apr2015/21apr2015
# Z Norgaard
#
# 21apr2015 change how insertions are represented in output.
# 10dec2015 made this version only identify indels
# 18mar2016 added multiple base changes back
#
#--> source('only_indels.R',echo=TRUE,max=Inf)
#
args = commandArgs(trailingOnly=T)

inmeta <- read.delim(args[1], header=T)
scientist <- args[2]

ff <- function(df, minv=1, maxv=5) df[minv:maxv, minv:maxv]
options(stringsAsFactors=FALSE) ### <<<<<<<<<<--------- NOTE!

library(ShortRead)
library(parallel)
library(VariantTools)
library(gmapR)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)
library(gtools)
library(stringr)

####################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # 
####################################################

# Today's date or what you will
if(FALSE) {
    dateTag = ''
} else {
    dateTag = tolower(format(Sys.Date(), "%d%b%Y"))
}
### END Required customization

inmeta <- inmeta[which(inmeta$Scientist == scientist),]

if(scientist=='Kevin') {
    file.dir = 'aligned_reads'
    reference = 'KH_corrected_15apr2016.fa'
    sp.start = 80
    sp.length = 2
    sp.plus = 49
}

if(scientist=='Olivier') {
    file.dir = 'aligned_reads'
    reference = 'OH_corrected_13jan2015.fa'
    sp.start = 120
    sp.length = 16
    sp.plus = 20
}

if(scientist=='Olivier_HbG') {
  file.dir = 'aligned_reads'
  reference = 'OH_HbG_16sep2016.fa'
  sp.start = 79
  sp.length = 20
  sp.plus = 10
}

if(scientist=='JS') {
    file.dir = 'sam_and_bam'
    reference = 'JS_background.fa'
    sp.start = 72
    sp.length = 6
    sp.plus = 1
}

if(!exists('file.dir')) stop('Customization failure!')

### Reference
refer = readFasta(reference)
al = ampliconLen = width(refer)

## mclapply cores
mycores = 5

### Two uppercase in a row function
  v = c('A', 'C', 'G', 'T')
  pe = permutations(4,2,v, repeats.allowed=TRUE)
  diplo = sapply(seq(nrow(pe)), function(k) paste(pe[k,],collapse=''))
  digrep <- function(y) any(sapply(diplo, function(XX) grepl(XX, y)))

### Filtered mutation frequencies
mm.freqFile = paste(scientist, '_filtered_spacer_indel_frequencies_',
                dateTag, '.txt', sep='')

hdr = c('Sample_name', 'Count', 'Mismatch_count', 'Match_frequency',
                                                  'Mismatch_frequency')
write(hdr, file=mm.freqFile, ncolumns=length(hdr), sep='\t')

### Stitched FASTQ files
bam = paste0(inmeta$File, '.bam')

for(abam in bam) {
  bamf = file.path(file.dir, abam)
  baseName = sub('.bam','', abam)

### Review stitched & filtered BAM file
  readGA = readGAlignments(bamf,
      param=ScanBamParam(what=c('seq'))) # qual and mapq removed for needle

  gadf = as.data.frame(readGA)
  gadf$seq = tolower(gadf$seq)
  fa = as.character(sread(refer))
  cg = gadf$cigar

### Start decomposing cigar
  cot = colSums(cigarOpTable(cigar(readGA)))
  cg.alpha = paste(names(cot[cot!=0]), sep='', collapse='') # 'MIDS'
  cgx = paste('[', cg.alpha, ']', sep='')

  break.cigar <- function(start, end) {
    cigar = cg[start:end]
#
    vals = strsplit(cg[start:end], cgx)
    cg.length = lapply(vals, as.integer)
#
    mids = strsplit(cg[start:end], '[[:digit:]]*')
    string = lapply(mids, function(x) x[seq(2, length(x), 2)])
#
    return(list(cigar=cigar, cg.length=cg.length, string=string))
  }

### Extract M, I, D, S cigar information (!!NB!! is that all the codes in KGH?)
  mids = break.cigar(1, length(readGA))

### Start constructing a cigar-split read set
  mids$kseq = mclapply(seq(length(mids$cigar)), function(k) {
    obj = lapply(mids, '[[', k)
    gs = gadf$seq[k]
#
    cs = cumsum(obj$cg.length)
    kstart = c(1, cs[-length(cs)]+1)
    kend = cs
    return(lapply(seq(length(obj$cg.length)),
        function(n) substr(gs, kstart[n], kend[n])))
  }, mc.cores=mycores)

### Re-construct the sequences
  recon = mclapply(seq(length(mids$cigar)), function(k) {
    obj = lapply(mids, '[[', k)
    iread = kread = ''  ### iread places inserts inline
    for(n in seq(length(obj$kseq))) {
        if(grepl('[SM]', obj$string[n])) {
            kread = paste(kread, obj$kseq[[n]], sep='')
            next
        }
        if(obj$string[n]=='I') {
            if(nchar(kread) >= sp.start - sp.plus &
               nchar(kread) <= sp.start + sp.length + sp.plus - 1) {
               iread = paste(iread, nchar(kread) + 1,
                   '.', toupper(obj$kseq[[n]]), ' ', sep='')
                len = nchar(kread)
                substr(kread,len,len) = toupper(substr(kread,len,len))
                if(n < length(obj$kseq)) substr(obj$kseq[[n+1]],1,1) =
                    toupper(substr(obj$kseq[[n+1]],1,1))
            }
            next
        }
        if(obj$string[n]=='D') {
            kread = paste(kread,
                paste(rep('-', obj$cg.length[[n]]), collapse=''),
                obj$kseq[[n]], sep='')
            next
        }
    }
    return(list(kread=kread, iread=iread)) 
  }, mc.cores=mycores)

  reads = unlist(sapply(recon, '[', 'kread'), use.names=FALSE)
  iread = unlist(sapply(recon, '[', 'iread'), use.names=FALSE)
  use <- nchar(reads) == nchar(fa)
  reads <- reads[use]
  iread <- iread[use]
  
### Make sure everything is the correct length
  stopifnot(all(nchar(reads) == nchar(fa)))
  
###
  fa.ok = tolower(fa)
  rf.ok = unlist(strsplit(fa.ok, ''))
  reads = unlist(mclapply(seq(length(reads)), function(k) {
    zed = unlist(strsplit(reads[k], ''))
    paste(ifelse(tolower(zed) != rf.ok & zed != '-', toupper(zed), zed), collapse='')
    }, mc.cores=mycores))

### Check for mutations in window
  subrf.ok = rf.ok[(sp.start - sp.plus) : (sp.start + sp.length + sp.plus - 1)]
  reads = unlist(mclapply(seq(length(reads)), function(k) {
    alph = unlist(strsplit(substring(reads[k], sp.start - sp.plus, 
                     sp.start + sp.length + sp.plus - 1), ''))
    # No mutations in window
    if (all(alph == subrf.ok)) {
      paste(fa.ok)
    # Mutation(s) in window
    } else {
      indel <- gregexpr('-+|[A-Z]+', reads[k])
      strs <- unlist(indel)
      # Verify at least one mutation is present
      stopifnot(all(strs != -1))
      if (length(strs) == 1) { # Only one mutation, logically must be in window
        paste(reads[k])
      } else { # Any number of mutations, logically at least one in window
        ends <- strs + attributes(indel[[1]])$match.length - 1
        chk <- logical()
        for (m in 1:length(strs)) {
          chk[m] <- (strs[m] > sp.start - sp.plus && 
                       strs[m] < sp.start + sp.length + sp.plus - 1) ||
            (ends[m] > sp.start - sp.plus &&
               ends[m] < sp.start + sp.length + sp.plus - 1) ||
            (sp.start > strs[m] && sp.start < ends[m])
        }
        sq <- ""
        st <- 1
        for (m in which(chk == TRUE)) {
          sq <- paste0(sq, substr(fa.ok, st, strs[m]-1), substr(reads[k], strs[m], ends[m]))
          st <- ends[m] + 1
        }
        sq <- paste0(sq, substr(fa.ok, st, nchar(fa.ok)))
        paste(sq)
      }
    }
  }, mc.cores=mycores))
  
### Organize mutations
  mu.all = data.frame(seq=names(table(reads)),
                count=as.vector(table(reads)))
  
  mu.all = mu.all[order(-mu.all$count),]
  rownames(mu.all) = NULL
  
### Compare highest frequency clone to reference
  #ok = mu.all$seq[1] # Should be FA, statistically speaking # 21apr2015 # 15jul2016 what's the point?, not technically true with high editing frequency
  #stopifnot(identical(ok, fa.ok))
  ok = fa.ok

  oks = unlist(strsplit(ok, ''))
  
  mu = mu.all[mu.all$count > 1,]
  mu <- mu[order(-mu$count),]
  rownames(mu) <- NULL
  
### Copy over insertions
  mu$iread <- ''
  two.mu = which(sapply(mu$seq, digrep, USE.NAMES = FALSE))
  for(k in two.mu) {
    mu$iread[k] = paste(unique(iread[reads==mu$seq[k]]), collapse='')
  }
    
### Classify Mutants
  insert = sapply(seq(nrow(mu)), function(k) # Insert counts
       length(unique(unlist(strsplit(gsub(' ', '', mu$iread[k]), '\\.[A-Z]*')))))
  
  # Added a better way to count deletions
  delete = sapply(seq(nrow(mu)), function(k)
    length(unlist(lapply(strsplit(mu$seq[k], split="[[:alpha:]]+"), 
                         function(x) {x[!x=='']}))))
  
  #Didn't just lower mu$seq in case of later compatibility issues
  mu$spacer_plus <- unlist(mclapply( seq(nrow(mu)), function(k) {
    zed = unlist(strsplit(mu$seq[k], ''))
    paste(ifelse(tolower(zed) == rf.ok, rf.ok, zed), collapse='')
  }, mc.cores=mycores))
  
  bpchange = sapply(seq(nrow(mu)), function(k)
    length(grep('[A-Z]', unlist(strsplit(mu$spacer_plus[k], ''))))
    )

  mu$changes = sapply(seq(nrow(mu)), function(k) {delete[k] + insert[k] + bpchange[k]})

  mu$spacer = substr(mu$spacer_plus, sp.start - sp.plus, sp.start + sp.length + sp.plus - 1)

### Indels or 2 changes in window
  mm.del <- grepl('-', mu$spacer_plus)
  mm.in <- mu$iread != ''
  mm.two <- mu$changes > 1
  mismatch = mm.del | mm.two | mm.in
  mu.2rule = rbind(mu[which(mu$changes == 0),], mu[mismatch,])
  ### "No mismatch" count includes mismatches eliminated by 2-mismatch rule
  mu.2rule$count[1] = mu.2rule$count[1] + sum(mu$count) - sum(mu.2rule$count)

### Ta-da
  mum = mu.2rule[, c('spacer_plus',  'count', 'spacer', 'changes', 'iread')]
  colnames(mum)[colnames(mum)=='iread'] = 'insert'
  #mum = mum[order(-mum$count),]

  mums = split(mum, mum$spacer)
  mumc = sapply(seq(length(mums)), function(k) sum(mums[[k]]$count))
  names(mumc) = names(mums)
  mumco = sort(mumc, decreasing=TRUE)
  mumcos = mums[names(mumco)]

  mumcoso = lapply(seq(length(mumcos)), function(k) {
    mumcos[[k]]$group = k
    return(mumcos[[k]])
  })

  mumcosorb = do.call('rbind', mumcoso)
  mumcosorb$freq = mumcosorb$count/sum(mumcosorb$count) # 21apr
  mumcosorb$freq = mumcosorb$count/sum(mu$count)
  mumcosorbo = mumcosorb[,
     c('spacer', 'spacer_plus','count','freq', 'changes','group', 'insert')]

  write.table(mumcosorbo, file=paste(baseName,
            '_spacer_plus', sp.plus, '_mutation_frequencies_',
            dateTag, '.txt', sep=''),
            quote=FALSE, sep='\t', row.names=FALSE)

### Summary info -- note details here!
#   o.total = sum(mumcosorbo$count)
#   o.mismatch = o.total - mumcosorbo$count[1]
#   o.mm.freq = o.mismatch/o.total
#   o.m.freq = 1 - o.mm.freq
#   out = c(baseName, o.total, o.mismatch, format(o.m.freq), format(o.mm.freq))
# 
#   #write(out, file=mm.freqFile, ncolumns=length(out), sep='\t', append=TRUE)
# 
# ### Now summarize by spacer
#   muml = split(mumcosorb, mumcosorb$spacer)
#   mumf = lapply(muml, function(df) 
#     data.frame(spacer=df$spacer[1], 
#                count=sum(df$count), freq=sum(df$count)/sum(mumco),
#                changes=df$changes[1],
#                group=df$group[1])
#     )
#   omum = do.call('rbind', mumf)
#   omum = omum[order(omum$group),]
#   rownames(omum) = NULL
# 
#   write.table(omum, file=paste(baseName,
#             '_spacer_mutation_frequencies_', dateTag, '.txt', sep=''),
#             quote=FALSE, sep='\t', row.names=FALSE)
#   print(baseName)
}
