#----------------------------------------------------------------------#
# Program: functions.R
# Author: Tim Simmons
# Date: 2019-02-14
# Purpose: Define general functions for Graham project
#----------------------------------------------------------------------#

# General utilities

first <- function(x) x[!is.na(x)][1]
avg <- function(x) mean(x, na.rm=TRUE)

as.set <- function(...) {
    x <- as.character(c(...))
    setNames(as.list(x), x)
}

# Data
read.tbl <- function(file, ...) {
    meta.file <- gsub(".csv$", "_labels.csv", file)

    if (!file.exists(meta.file)) {
        tbl <- read.csv(
            file,
            stringsAsFactors=FALSE,
            na.strings=c("", "Not Available", "NotAvailable", "NA"))
        variables <- names(tbl)
        types <- sapply(tbl, typeof)
        types[types == "double"] <- "numeric"
        substring(types, 1, 1) <- toupper(substring(types, 1, 1))

        meta <- data.frame(
            Variable=variables,
            Types=types,
            Label="",
            stringsAsFactors=FALSE)
        write.csv(meta, file=meta.file, row.names=FALSE)
    } else {
        meta <- read.csv(meta.file, stringsAsFactors=FALSE)
    }

    colClasses <- tolower(meta$Type)
    names(colClasses) <- meta$Variable

    irf <- read.csv(
        file,
        stringsAsFactors=FALSE,
        na.strings=c("", "Not Available", "NotAvailable", "NA"),
        colClasses=colClasses,
        ...
    )
    irf
}

read.main <- function() {
    main <- read.tbl("../data/irfmain.csv")

    main$region <- factor(
        main$region,
        as.character(1:10),
        c("Boston", "New York"   , "Philadelphia", "Atlanta", "Chicago",
          "Dallas", "Kansas City", "Denver", "San Francisco", "Seattle"))

    main$census <- factor(
        main$census,
        1:9,
        c("New England", "Middle Atlantic", "East North Central",
          "West North Central", "South Atlantic", "East South Central",
          "West South Central", "Mountain", "Pacific"))

    main$rural <- factor(
        main$rural,
        0:1,
        c("Urban", "Rural"))

    main$ownership <- factor(
        main$ownership,
        c("G", "N", "P"),
        c("Government", "Non-profit", "Profit"))

    main$teaching.yn <- factor(
        main$teaching.yn,
        0:1,
        c("No", "Yes"))

    main$freestanding <- factor(
        main$freestanding,
        0:1,
        c("No", "Yes"))

    main
}


is.good.grade <- function(x) {
    y <- integer(length(x))
    y[] <- NA_integer_

    y[grepl("^Better", x)] <- 1L
    y[grepl("^No Different|^Worse", x)] <- 0L
    y
}

std <- function(df) do.call(cbind, lapply(df, function(x) (x-mean(x))/sd(x)))

format.p.value <- function(p) {
    s <- rep("NULL", length(p))
    s[p >= 0.01] <- sprintf("%.2f", p[which(p >= 0.01)])
    s[p >= 0.001 & p < 0.01] <- sprintf("%.3f", p[which(p >= 0.001 & p < 0.01)])
    s[p < 0.001] <- "< 0.001"
    s
}

# Generic code to estimate mean responses on new data from either lm or
# lme objects
augur <- function(obj, data) {

    fml <- as.formula(paste(as.character(formula(obj))[-2], collapse=""))
    X <- model.matrix(fml, model.frame(fml, data))

    b <- X %*% coef(summary(obj))[, 1]
    V <- X %*% vcov(obj) %*% t(X)

    list(fit=b, se.fit=sqrt(diag(V)), vcov=V)
}

augur.compare <- function(obj, alpha=0.05, method="bonferroni") {
    b <- obj$fit
    p <- length(b)
    if (method == "bonferroni") alpha <- alpha*2/p/(p-1)
    z.alpha <- qnorm(1 - alpha/2)

    a <- array(0, dim=c(1L, p))
    z <- array(2, dim=c(p, p))

    for (i in 1L:(p-1L)) {
        for (j in (i+1L):p) {
            a[] <- 0
            a[, i] <- 1
            a[, j] <- -1

            x <- a %*% b
            se <- sqrt(a %*% obj$vcov %*% t(a))
            z[i, j] <- z[j, i] <- 2*pnorm(abs(x)/se, lower.tail=FALSE)
        }
    }

    list(p.values=z, b.values=z*(p-1)*p/2, signficance=(z < alpha))
}

profile <- function(df) {
    col <- names(df)
    summary <- lapply(col, function(x) data.frame(
        column=x,
        n=sum(!is.na(df[[x]])),
        n.missing=sum(is.na(df[[x]])),
        n.unique=length(unique(df[[x]]))
    ))
    do.call(rbind, summary)
}



#----------------------------------------------------------------------#
# Closure: cld
# Source: multcomp_1.4-10.tar.gz
# Notes:
#   Functions from the multcomp package appear unchanged, but wrapped in
# a closure to keep the workspace clean.
#----------------------------------------------------------------------#

cld <- (function() {
# Function implements the insert-absorb (sweep) heuristic of Piepho 2004:
# "An Algorithm for a Letter-Based Representation of All-Pairwise Comparisons"
#
# x         ... vector of logicals indicating significant comparisons with hyphenated
#               names e.g. A-B, treatmentA-treatmentB, ...
# Letters   ... a set of user defined letters { default is Letters=c(letters, LETTERS) }
# separator ... a separating character used to produce a sufficiently large set of
#               characters for a compact letter display (default is separator=".") in case
#               the number of letters required exceeds the number of letters available
# Decreasing ... Inverse the order of the letters

insert_absorb <- function( x, Letters=c(letters, LETTERS), separator=".", decreasing = FALSE,
                           comps = NULL, lvl_order){

  obj_x <- deparse(substitute(x))
  if (is.null(comps)) {
      namx <- names(x)
      namx <- gsub(" ", "", names(x))
      if(length(namx) != length(x))
          stop("Names required for ", obj_x)
      split_names <- strsplit(namx, "-")
      stopifnot( sapply(split_names, length) == 2 )
      comps <- t(as.matrix(as.data.frame(split_names)))
  }
  rownames(comps) <- names(x)
  lvls <- lvl_order
  n <- length(lvls)
  lmat <- array(TRUE, dim=c(n,1), dimnames=list(lvls, NULL) )

  if( sum(x) == 0 ){                                                        # no differences
    ltrs <- rep(get_letters(1, Letters=Letters, separator=separator), length(lvls) )
    names(ltrs) <- lvls
    colnames(lmat) <- ltrs[1]
    msl <- ltrs
    ret <- list(Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat)
    class(ret) <- "multcompLetters"
    return(ret)
  }
  else{
    signifs <- comps[x,,drop=FALSE]

    absorb <- function(m){
      for(j in 1:(ncol(m)-1)){
        for(k in (j+1):ncol(m)){
          if( all(m[which(m[,k]),k] & m[which(m[,k]),j]) ){                 # column k fully contained in column j
            m <- m[,-k, drop=FALSE]
            return(absorb(m))
          }
          else if( all(m[which(m[,j]),k] & m[which(m[,j]),j]) ){            # column j fully contained in column k
            m <- m[,-j, drop=FALSE]
            return(absorb(m))
          }
        }
      }
      return(m)
    }
    for( i in seq_len(nrow(signifs)) ){                                            # insert
      tmpcomp <- signifs[i,]
      wassert <- which(lmat[tmpcomp[1],] & lmat[tmpcomp[2],])               # which columns wrongly assert nonsignificance
      if(any(wassert)){
        tmpcols <- lmat[,wassert,drop=FALSE]
        tmpcols[tmpcomp[2],] <- FALSE
        lmat[tmpcomp[1],wassert] <- FALSE
        lmat <- cbind(lmat, tmpcols)
        colnames(lmat) <- get_letters( ncol(lmat), Letters=Letters,
                                       separator=separator)
        if(ncol(lmat) > 1){                                                 # absorb columns if possible
          lmat <- absorb(lmat)
          colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                         separator=separator )
        }
      }
    }
  }
  lmat <- lmat[,order(apply(lmat, 2, sum))]
  lmat <- sweepLetters(lmat)                                                                  # 1st
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x))))))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  lmat <- lmat[,order(apply(lmat, 2, sum))]                                                   # 2nd sweep
  lmat <- sweepLetters(lmat)
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x)))),
                           decreasing = decreasing))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  ltrs <- apply(lmat,1,function(x) return(paste(names(x)[which(x)], sep="", collapse="") ) )
  msl <- matrix(ncol=ncol(lmat), nrow=nrow(lmat))                                             # prepare monospaced letters
  for( i in 1:nrow(lmat) ){
    msl[i,which(lmat[i,])] <- colnames(lmat)[which(lmat[i,])]
    absent <- which(!lmat[i,])
    if( length(absent) < 2 ){
      if( length(absent) == 0 )
        next
      else{
        msl[i,absent] <- paste( rep(" ", nchar(colnames(lmat)[absent])), collapse="" )
      }
    }
    else{
      msl[i,absent] <- unlist( lapply( sapply( nchar(colnames(lmat)[absent]),
                                               function(x) return(rep( " ",x)) ),
                                       paste, collapse="") )
    }
  }
  msl <- apply(msl, 1, paste, collapse="")
  names(msl) <- rownames(lmat)
  ret <- list( Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat,
               aLetters = Letters, aseparator = separator )
  class(ret) <- "multcompLetters"
  return(ret)
}


# All redundant letters are swept out without altering the information within a LetterMatrix.
#
# mat         ... a LetterMatrix as produced by function insert_absorb()
# start.col   ... either a single integer specifying the column to start with or a vector
#                 of max. length equal to ncol(mat) specifying the column order to be used.
# Letters     ... a set of user defined letters { default is Letters=c(letters, LETTERS) }
# separator   ... a separating character used to produce a sufficiently large set of
#                 characters for a compact letter display (default is separator=".") in case
#                 the number of letters required exceeds the number of letters available

sweepLetters <- function(mat, start.col=1, Letters=c(letters, LETTERS), separator="."){

  stopifnot( all(start.col %in% 1:ncol(mat)) )
  locked <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))          # 1 indicates that another letter dependes on this entry
  cols <- 1:ncol(mat)
  cols <- cols[c( start.col, cols[-start.col] )]
  if( any(is.na(cols) ) )
    cols <- cols[-which(is.na(cols))]

  for( i in cols){
    tmp <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))
    tmp[which(mat[,i]),] <- mat[which(mat[,i]),]                        # get items of those rows which are TRUE in col "i"
    one <- which(tmp[,i]==1)

    if( all(apply(tmp[,-i,drop=FALSE], 1, function(x) return( any(x==1) ))) ){     # there is at least one row "l" where mat[l,i] is the only item which is TRUE i.e. no item can be removed in this column
      next
    }
    for( j in one ){                                                    # over all 1's
      if( locked[j,i] == 1 ){                                           # item is locked
        next
      }
      chck <- 0
      lck <- list()
      for( k in one ){
        if( j==k ){
          next
        }
        else{                                                           # pair j-k
          rows <- tmp[c(j,k),]
          dbl <- rows[1,] & rows[2,]
          hit <- which(dbl)
          hit <- hit[-which(hit==i)]
          dbl <- rows[1,-i,drop=FALSE] & rows[2,-i,drop=FALSE]
          if( any(dbl) ){
            chck <- chck + 1
            lck[[chck]] <- list(c(j,hit[length(hit)]), c(k,hit[length(hit)]))      # record items which have to be locked, use last column if multiple hits
          }
        }
      }
      if( (chck == (length(one)-1)) && chck != 0 ){                     # item is redundant
        for( k in 1:length(lck) ){                                      # lock items
          locked[ lck[[k]][[1]][1], lck[[k]][[1]][2] ] <- 1
          locked[ lck[[k]][[2]][1], lck[[k]][[2]][2] ] <- 1
        }
        mat[j,i] <- FALSE                                               # delete redundant entry
      }
    }
    if(all(mat[,i]==FALSE)){                                           # delete column where each entry is FALSE and restart
      mat <- mat[,-i,drop=FALSE]
      colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
      return(sweepLetters(mat, Letters=Letters, separator=separator))
    }
  }
  onlyF <- apply(mat, 2, function(x) return(all(!x)))
  if( any(onlyF) ){                                                     # There are columns with just FALSE entries
    mat <- mat[,-which(onlyF),drop=FALSE]
    colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
  }
  return( mat )
}

# Create a set of letters for a letter display. If "n" exceeds the number of letters
# specified in "Letters", they are recycled with one or more separating character(s)
# preceding each recycled letter.
# e.g. get_letters(10, Letters=letters[1:4]) produces:  "a"   "b"   "c"   "d"   ".a"  ".b"  ".c"  ".d"  "..a" "..b"
#
# n             ... number of letters
# Letters       ... the set of characters to be used
# separator     ... a character to be used as separator e.g.
#                   n=5, Letters=c("a","b") => "a", "b", ".a", ".b", "..a"

get_letters <- function( n, Letters=c(letters, LETTERS), separator="." ){

  n.complete <- floor(n / length(Letters))        # number of complete sets of Letters
  n.partial <- n %% length(Letters)               # number of additional Letters
  lett <- character()
  separ=""
  if( n.complete > 0 ){
    for( i in 1:n.complete ){
      lett <- c(lett, paste(separ, Letters, sep="") )
      separ <- paste( separ, separator, sep="" )
    }
  }
  if(n.partial > 0 )
    lett <- c(lett, paste(separ, Letters[1:n.partial], sep="") )
  return(lett)
}

list(insert_absorb=insert_absorb)
})()