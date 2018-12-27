# Playing around with Katz Back off algorithm example 
# from https://rpubs.com/mszczepaniak/predictkbo3model
setwd('~/GitRepos/InitialCapstone/')
library(data.table)
library(stringr)
ltcorpus<-readLines("thelittlecorpus.txts")

## Returns a named vector of n-grams and their associated frequencies
## extracted from the character vector dat.
##
## ng - Defines the type of n-gram to be extracted: unigram if ng=1,
##      bigram if ng=2, trigram if n=3, etc.
## dat - Character vector from which we want to get n-gram counts.
## ignores - Character vector of words (features) to ignore from frequency table
## sort.by.ngram - sorts the return vector by the names
## sort.by.freq - sorts the return vector by frequency/count
getNgramFreqs <- function(ng, dat, ignores=NULL,
                          sort.by.ngram=TRUE, sort.by.freq=FALSE) {
  # http://stackoverflow.com/questions/36629329/
  # how-do-i-keep-intra-word-periods-in-unigrams-r-quanteda
  if(is.null(ignores)) {
    dat.dfm <- dfm(dat, ngrams=ng, toLower = FALSE, removePunct = FALSE,
                   what = "fasterword", verbose = FALSE)
  } else {
    dat.dfm <- dfm(dat, ngrams=ng, toLower = FALSE, ignoredFeatures=ignores,
                   removePunct = FALSE, what = "fasterword", verbose = FALSE)
  }
  rm(dat)
  # quanteda docfreq will get the document frequency of terms in the dfm
  ngram.freq <- docfreq(dat.dfm)
  if(sort.by.freq) { ngram.freq <- sort(ngram.freq, decreasing=TRUE) }
  if(sort.by.ngram) { ngram.freq <- ngram.freq[sort(names(ngram.freq))] }
  rm(dat.dfm)
  
  return(ngram.freq)
}

## Returns a 2 column data.table. The first column: ngram, contains all the
## unigrams, bigrams, or trigrams in the corpus depending on whether
## ng = 1, 2, or 3 respectively. The second column: freq, contains the
## frequency or count of the ngram found in linesCorpus.
##
## ng - Defines the type of n-gram to be extracted: unigram if ng=1,
##      bigram if ng=2, trigram if n=3, etc.
## linesCorpus - character vector: each element is a line from a corpus file
## prefixFilter - character vector: If not NULL, tells the function to return
##                only rows where the ngram column starts with prefixFilter.
##                If NULL, returns all the ngram and count rows.
getNgramTables <- function(ng, linesCorpus, prefixFilter=NULL) {
  ngrams <- getNgramFreqs(ng, linesCorpus)
  ngrams_dt <- data.table(ngram=names(ngrams), freq=ngrams)
  if(length(grep('^SOS', ngrams_dt$ngram)) > 0) {
    ngrams_dt <- ngrams_dt[-grep('^SOS', ngrams_dt$ngram),]
  }
  if(!is.null(prefixFilter)) {
    regex <- sprintf('%s%s', '^', prefixFilter)
    ngrams_dt <- ngrams_dt[grep(regex, ngrams_dt$ngram),]
  }
  
  return(ngrams_dt)
}

unigs <- getNgramTables(1, ltcorpus)
bigrs <- getNgramTables(2, ltcorpus)
trigs <- getNgramTables(3, ltcorpus)

gamma2 <- 0.5  # bigram discount
gamma3 <- 0.5  # trigram discount
bigPre <- 'sell_the'

## Returns a two column data.frame of observed trigrams that start with the
## bigram prefix (bigPre) in the first column named ngram and
## frequencies/counts in the second column named freq. If no observed trigrams
## that start with bigPre exist, an empty data.frame is returned.
##
## bigPre -  single-element char array of the form w2_w1 which are the first 
##           two words of the trigram we are predicting the tail word of
## trigrams - 2 column data.frame or data.table. The first column: ngram,
##            contains all the trigrams in the corpus. The second column:
##            freq, contains the frequency/count of each trigram.
getObsTrigs <- function(bigPre, trigrams) {
  trigs.winA <- data.frame(ngrams=vector(mode = 'character', length = 0),
                           freq=vector(mode = 'integer', length = 0))
  regex <- sprintf("%s%s%s", "^", bigPre, "_")
  trigram_indices <- grep(regex, trigrams$ngram)
  if(length(trigram_indices) > 0) {
    trigs.winA <- trigrams[trigram_indices, ]
  }
  
  return(trigs.winA)
}

## Returns a two column data.frame of observed trigrams that start with bigram
## prefix bigPre in the first column named ngram and the probabilities
## q_bo(w_i | w_i-2, w_i-1) in the second column named prob calculated from
## eqn 12. If no observed trigrams starting with bigPre exist, NULL is returned.
##
## obsTrigs - 2 column data.frame or data.table. The first column: ngram,
##            contains all the observed trigrams that start with the bigram
##            prefix bigPre which we are attempting to the predict the next
##            word of in a give phrase. The second column: freq, contains the
##            frequency/count of each trigram.
## bigrs - 2 column data.frame or data.table. The first column: ngram,
##         contains all the bigrams in the corpus. The second column:
##         freq, contains the frequency/count of each bigram.
## bigPre -  single-element char array of the form w2_w1 which are first two
##           words of the trigram we are predicting the tail word of
## triDisc - amount to discount observed trigrams
getObsTriProbs <- function(obsTrigs, bigrs, bigPre, triDisc=0.5) {
  if(nrow(obsTrigs) < 1) return(NULL)
  obsCount <- filter(bigrs, ngram==bigPre)$freq[1]
  obsTrigProbs <- mutate(obsTrigs, freq=((freq - triDisc) / obsCount))
  colnames(obsTrigProbs) <- c("ngram", "prob")
  
  return(obsTrigProbs)
}

obs_trigs <- getObsTrigs(bigPre, trigs)  # get trigrams and counts
# convert counts to probabilities
qbo_obs_trigrams <- getObsTriProbs(obs_trigs, bigrs, bigPre, gamma3)
qbo_obs_trigrams

## Returns a character vector which are the tail words of unobserved trigrams
## that start with the first two words of obsTrigs (aka the bigram prefix).
## These are the words w in the set B(w_i-2, w_i-1) as defined in the section
## describing the details of equation 17.
##
## obsTrigs - character vector of observed trigrams delimited by _ of the form:
##            w3_w2_w1 where w3_w2 is the bigram prefix
## unigs - 2 column data.frame of all the unigrams in the corpus:
##         ngram = unigram
##         freq = frequency/count of each unigram
getUnobsTrigTails <- function(obsTrigs, unigs) {
  obs_trig_tails <- str_split_fixed(obsTrigs, "_", 3)[, 3]
  unobs_trig_tails <- unigs[!(unigs$ngram %in% obs_trig_tails), ]$ngram
  return(unobs_trig_tails)
}

unobs_trig_tails <- getUnobsTrigTails(obs_trigs$ngram, unigs)
unobs_trig_tails

## Returns the total probability mass discounted from all observed bigrams
## calculated from equation 14.  This is the amount of probability mass which
## is redistributed to UNOBSERVED bigrams. If no bigrams starting with
## unigram$ngram[1] exist, 0 is returned.
##
## unigram - single row, 2 column frequency table. The first column: ngram,
##           contains the w_i-1 unigram (2nd word of the bigram prefix). The
##           second column: freq, contains the frequency/count of this unigram.
## bigrams - 2 column data.frame or data.table. The first column: ngram,
##           contains all the bigrams in the corpus. The second column:
##           freq, contains the frequency or count of each bigram.
## bigDisc - amount to discount observed bigrams
getAlphaBigram <- function(unigram, bigrams, bigDisc=0.5) {
  # get all bigrams that start with unigram
  regex <- sprintf("%s%s%s", "^", unigram$ngram[1], "_")
  bigsThatStartWithUnig <- bigrams[grep(regex, bigrams$ngram),]
  if(nrow(bigsThatStartWithUnig) < 1) return(0)
  alphaBi <- 1 - (sum(bigsThatStartWithUnig$freq - bigDisc) / unigram$freq)
  
  return(alphaBi)
}

unig <- str_split(bigPre, "_")[[1]][2]
unig <- unigs[unigs$ngram == unig,]
alpha_big <- getAlphaBigram(unig, bigrs, gamma2)
alpha_big