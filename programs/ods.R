#----------------------------------------------------------------------#
# Program: ods.R
# Author: Tim Simmons
# Date: 2018-10-19
# Purpose: Defines functions for putting output to an HTML document
# Depends: base64enc
# TODO:
#  1. Functions for tabulation
#  2. Proper, consistent escaping
#  3. Copy css and js source files into HTML document
#----------------------------------------------------------------------#


html <- function(htmlfile) {

    doc <- file(htmlfile, "w")

    `%<<%` <- function(con, x) write(x, con, append=TRUE)

    cssfile <- gsub("\\.html$", ".css", basename(htmlfile))
    jsfile <- gsub("\\.html$", ".js", basename(htmlfile))

    doc %<<% "<!DOCTYPE html>\n"
    doc %<<% "<html lang=\"en\">"
    doc %<<% "<head>\n<meta charset=\"utf-8\">\n"
    doc %<<% sprintf("<link rel=\"stylesheet\" href=\"%s\">\n", cssfile)
    doc %<<% "</head>\n<body>\n"

    nfig <- 0L
    ntab <- 0L
    caption <- ""

    title <- function(s) caption <<- s

    escape <- function(text) text

    tag <- function(t) {
        function(s) doc %<<% sprintf("<%s>%s</%s>\n", t, escape(s), t)
    }

    put <- function(s) doc %<<% s

    as.string <- function(...) {
        con <- textConnection("s", open="w", local=TRUE)
        sink(con)
        base::print(...)
        sink()
        base::close(con)
        paste(s, collapse="\n")
    }

    # Print object as it would appear in the R console
    raw <- function(...) tag("pre")(as.string(...))

    # Print an html representation of the object (if possible)
    print <- function(x, ...) {
        cls <- class(x)
        doc %<<% if (!is.null(cls) && cls != "") sprintf("<section class=\"%s\">", gsub("[.]", "-", cls)) else "<section>"

        tabulator.classes <- c(
            "tabulator.n",
            "tabulator.table",
            "tabulator.mean",
            "tabulator.median")

        if (is.data.frame(x)) .print.data.frame(x, ...)
        else if (class(x) == "matrix") .print.matrix(x)
        else if (class(x) == "tabulator") .print.tabulator(x)
        else if (class(x) %in% tabulator.classes) .print.matrix(unclass(x))
        else raw(x)

        doc %<<% "</section>"
    }

    .print.matrix <- function(x) .print.data.frame(as.data.frame(x))
    .print.data.frame <- function(x, row.names=TRUE) {
        doc %<<% "<table>"


        if (caption != "") {
            ntab <<- ntab + 1L
            text <- sprintf("Table %d: %s", ntab, caption)
            tag("caption")(text)
            caption <<- ""
        }

        doc %<<% "<thead>"
        doc %<<% "<tr>"
        if (row.names) doc %<<% "<td></td>"
        doc %<<% paste(sprintf("<th>%s</th>", colnames(x)), collapse="")
        doc %<<% "</tr>"
        doc %<<% "</thead>"


        doc %<<% "<tbody>"
        if (row.names) row <- sprintf("<th>%s</th>", rownames(x))
        data <- do.call(
            paste0,
            lapply(x, function(x) {
                type <- if (class(x) %in% c("integer", "numeric")) "n" else "c"
                sprintf("<td class=\"td-r-%s\">%s</td>", type, as.character(x))
        }))
        n <- nrow(x)
        for (i in 1:n) {
            doc %<<% "<tr>"
            if (row.names) doc %<<% row[i]
            doc %<<% data[i]
            doc %<<% "</tr>"
        }
        doc %<<% "</tbody>"

        doc %<<% "</table>"
    }
    .print.tabulator <- function(x) lapply(x$result(), print)

    # Copy plot in current device into the html document
    png <- function(...) {
        tmp <- tempfile()

        dev.copy(function() grDevices::png(tmp, ...))
        dev.off()

        img64 <- base64enc::base64encode(tmp)

        doc %<<% "<figure>"
        doc %<<% sprintf("<img src=\"data:image/png;base64, %s\"/>", img64)

        if (caption != "") {
            nfig <<- nfig + 1L
            text <- sprintf("Figure %d: %s", nfig, caption)
            tag("figcaption")(text)
            caption <<- ""
        }

        doc %<<% "</figure>"

        unlink(tmp)
    }

    # Close the html document
    close <- function() {
        doc %<<% sprintf("<script src=\"%s\"></script>\n", jsfile)
        doc %<<% "</body></html>"
        base::close(doc)
    }

    list(
        h1 = tag("h1"),
        h2 = tag("h2"),
        h3 = tag("h3"),
        h4 = tag("h4"),
        h5 = tag("h5"),
        h6 = tag("h6"),
        title = title,
        p = tag("p"),
        put = put,
        raw = raw,
        print = print,
        png = png,
        close = close
    )
}


tabulator <- function(data, across) {

    term.variables <- all.names(attr(terms(across), "variables"))[-1]
    for (t in term.variables) {
        if (!is.factor(data[[t]])) data[[t]] <- factor(data[[t]])
        levels(data[[t]]) <- paste0("/", levels(data[[t]]))
    }

    contr <- function(x) contr.treatment(levels(x), contrasts=FALSE)
    contr.arg <- lapply(data[term.variables], contr)

    x <- model.matrix(across, data=data, contrasts.arg=contr.arg)
    rows <- rownames(x)

    result <- list()
    i <- 1L

    n <- function() {
        counts <- colSums(x)
        dim(counts) <- c(1L, length(counts))
        colnames(counts) <- colnames(x)
        class(counts) <- "tabulator.n"
        result[[i]] <<- counts
        i <<- i + 1L
        counts
    }

    table <- function(var, n="%d", pct="%.1f%%") {
        f <- data[rows, var]
        groups <- split(rows, f)
        counts <- lapply(groups, function(g) colSums(x[g, , drop=FALSE]))
        counts <- do.call(rbind, counts)
        rownames(counts) <- paste(var, rownames(counts), sep="/")

        spct <- sprintf(pct, counts/rep(colSums(counts), each=nrow(counts))*100)
        dim(spct) <- dim(counts)
        rownames(spct) <- rownames(counts)

        counts <- rbind(counts, spct)

        class(counts) <- "tabulator.table"
        result[[i]] <<- counts
        i <<- i + 1L
        counts
    }

    mean <- function(var, n="%d", mean="%.1f", sd="%.1f") {
        fn <- n
        fmean <- mean
        fsd <- sd

        v <- data[rows, var]
        t <- do.call(cbind, lapply(colnames(x), function(n) {
            y <- v[as.logical(x[, n])]
            c("n"=sprintf(fn, sum(!is.na(y))),
              "mean"=sprintf(fmean, base::mean(y, na.rm=TRUE)),
              "sd"=sprintf(fsd, stats::sd(y, na.rm=TRUE)))
        }))
        rownames(t) <- paste(var, rownames(t), sep="/")
        colnames(t) <- colnames(x)
        class(t) <- "tabulator.mean"
        result[[i]] <<- t
        i <<- i + 1L
        t

    }

    median <- function(var, median="%.1f", q1="%.1f", q3="%.1f") {
        fmedian <- median
        fq1 <- q1
        fq3 <- q3

        v <- data[rows, var]
        t <- do.call(cbind, lapply(colnames(x), function(n) {
            y <- v[as.logical(x[, n])]
            s <- summary(y)
            c("median"=sprintf(fmedian, s["Median"]),
              "q1"=sprintf(fq1, s["1st Qu."]),
              "q3"=sprintf(fq3, s["3rd Qu."]))
        }))
        rownames(t) <- paste(var, rownames(t), sep="/")
        colnames(t) <- colnames(x)
        class(t) <- "tabulator.median"
        result[[i]] <<- t
        i <<- i + 1L
        t

    }

    structure(list(
        x=x,
        n=n,
        table=table,
        mean=mean,
        median=median,
        result=function() return(result)
    ), class="tabulator")
}