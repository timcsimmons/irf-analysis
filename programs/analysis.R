#----------------------------------------------------------------------#
# Program: analysis.R
# Author: Tim Simmons
# Date: 2019-02-14
# Purpose: Perform analyses for Graham project
#----------------------------------------------------------------------#

source("functions.R")
source("ods.R")

library(nlme)


output <- html("../output/output.html")

main <- read.main()
main$n.discharge100 <- main$n.discharge/100


facility.chrs <- c(
    "n.discharge",
    "cost.charge.ratio",
    "wage.index",
    "rural",
    "dsh",
    "census",
    "ownership",
    "teaching.yn",
    "freestanding",
    "population.severity")

process.measures <- c(
    "pt.flu",
    "functional",
    "hc.flu",
    "mspb")

outcome.measures <- c(
    "bedsores.v1",
    "bedsores.v2",
    "cauti",
    "readmit30",
    "falls",
    "cdi",
    "home")

km.outcome.measures <- c(
    "cauti",
    "readmit30",
    "cdi",
    "home")

outcome.grades <- c(
    "cauti.good",
    "readmit30.good",
    "cdi.good",
    "home.good")


#----------------------------------------------------------------------#
# List of process and outcome measures
#----------------------------------------------------------------------#

output$h1("Performance measures")

output$title("CMS codes and descriptions of the performance measures and years when they are available")

irfmeasures <- read.tbl("../data/irfmeasures.csv")

irfmeasures$years[is.na(irfmeasures$years)] <- ""
irfmeasures[["y2016"]] <- ifelse(grepl("2016", irfmeasures$year), "*", "")
irfmeasures[["y2017"]] <- ifelse(grepl("2017", irfmeasures$year), "*", "")
irfmeasures[["y2018"]] <- ifelse(grepl("2018", irfmeasures$year), "*", "")

irfmeasures <- do.call(rbind,
    lapply(split(irfmeasures, irfmeasures$type), function(x) {
        rbind(
            data.frame(cms.id="",
                       description=x$type[1],
                       y2016="",
                       y2017="",
                       y2018="",
                       stringsAsFactors=FALSE),
            x[c("cms.id", "description", "y2016", "y2017", "y2018")]
        )
    })[c(3, 2, 1)]
)

row.names(irfmeasures) <- NULL
output$print(irfmeasures, row.names=FALSE)


#----------------------------------------------------------------------#
# Descriptive statistics
#----------------------------------------------------------------------#

output$h1("Descriptive statistics")

output$title("Descriptive summaries of facility characteristics and performance measures by year")
t0 <- tabulator(main, ~ -1 + year)
t0$n()
t0$median("n.discharge", median="%.0f", q1="%.0f", q3="%.0f")
t0$mean("cost.charge.ratio", mean="%.2f", sd="%.2f")
t0$mean("wage.index", mean="%.2f", sd="%.2f")
t0$table("rural")
t0$median("dsh")
t0$table("census")
t0$table("ownership")
t0$table("teaching.yn")
t0$table("freestanding")
t0$mean("population.severity", mean="%.2f", sd="%.2f")

t0$mean("bedsores.v1", mean="%.2f", sd="%.2f")
t0$mean("bedsores.v2", mean="%.2f", sd="%.2f")
t0$mean("pt.flu"     , mean="%.2f", sd="%.2f")
t0$mean("cauti"      , mean="%.2f", sd="%.2f")
t0$mean("readmit30"  , mean="%.2f", sd="%.2f")
t0$mean("functional" , mean="%.2f", sd="%.2f")
t0$mean("falls"      , mean="%.2f", sd="%.2f")
t0$mean("cdi"        , mean="%.2f", sd="%.2f")
t0$mean("hc.flu"     , mean="%.2f", sd="%.2f")
t0$mean("home"       , mean="%.2f", sd="%.2f")
t0$mean("mspb"       , mean="%.2f", sd="%.2f")

t0$table("cauti.grade")
t0$table("readmit30.grade")
t0$table("cdi.grade")
t0$table("home.grade")

output$print(t0)


#----------------------------------------------------------------------#
# Aim 1: Association between process and outcome measures
#----------------------------------------------------------------------#


# Spearman correlation for each year (esp. 2018)
t1 <- lapply(split(main, main$year), function(u) {

    r <- array("", dim=c(length(outcome.measures),
                         length(process.measures)))

    for (i in seq_along(outcome.measures)) {
        for (j in seq_along(process.measures)) {

            y <- u[[outcome.measures[i]]]
            x <- u[[process.measures[j]]]
            if (sum(!is.na(x) & !is.na(y)) <= 2) next

            r0 <- cor.test(x, y, method="spearman", exact=FALSE)
            r[i, j] <- sprintf(
                "%.2f (%s)",
                r0$estimate, format.p.value(r0$p.value))
        }
    }

    t <- cbind(as.character(u$year[1]), outcome.measures, r)
    colnames(t) <- c("year", "outcome.measure", process.measures)
    t
})

output$h1("Correlation between process and outcome measures")
output$h2("Pairwise associations")
output$title("Spearman rank correlations (p-values) between process and outcome measures")

t1 <- do.call(cbind, t1[c("2017", "2018")])
colnames(t1) <- c(
    "year",
    "outcome.measure",
    paste0(process.measures, ":2017"),
    "year:2018",
    "outcome.measure:2018",
    paste0(process.measures, ":2018"))
t1 <- t1[, c("outcome.measure",
             paste0(rep(process.measures, each=2), c(":2017", ":2018")))]

output$print(t1)

# Canonical correlation (esp. 2018)
x <- as.matrix(main[main$year == 2018L, process.measures])[, -c(1, 3)]
y <- as.matrix(main[main$year == 2018L, outcome.measures])[, -c(1, 4)]
i <- complete.cases(x) & complete.cases(y)
rho <- cancor(x[i, ], y[i, ])
rho.p <- CCP::p.asym(rho$cor, sum(i), ncol(x), ncol(y), tstat="Wilks")

rx <- do.call(cbind, lapply(1:ncol(x), function(u) rank(x[i, u])))
ry <- do.call(cbind, lapply(1:ncol(y), function(u) rank(y[i, u])))
rrho <- cancor(rx, ry)
rrho.p <- CCP::p.asym(rrho$cor, sum(i), ncol(x), ncol(y), tstat="Wilks")

output$h2("Canonical correlations")
output$put(sprintf("<p>Canonical correlation analysis was applied to the %d facilities with complete process and outcome measures in 2018. For the unranked measures, the largest canonical correlation was %.2f (p %s by Wilk&apos;s lambda). For the ranked measures, the largest correlation was %.2f (p %s by Wilk&apos;s lambda).</p>\n", sum(i), rho$cor[1], format.p.value(rho.p$p.value[1]), rrho$cor[1], format.p.value(rrho.p$p.value[1])))


#----------------------------------------------------------------------#
# Aim 2: Compare high and low performing facilities
#----------------------------------------------------------------------#

output$h1("Associations between outcome grades and facility characteristics")
output$h2("k means clustering")


# k-means on measures with grades (complete cases) [2 clusters]
kmain <- do.call(
    rbind,
    lapply(split(main, main$facility.id), function(x) {
        data.frame(
            facility.id   = x$facility.id[1],
            cauti.avg     = avg(x$cauti),
            readmit30.avg = avg(x$readmit30),
            cdi.avg       = avg(x$cdi),
            home.avg      = avg(x$home),

            cost.charge.ratio.avg = avg(x$cost.charge.ratio),
            wage.index.avg= avg(x$wage.index),
            dsh.avg       = avg(x$dsh),
            population.severity.avg = avg(x$population.severity),

            rural.1       = first(x$rural),
            census.1      = first(x$census),
            ownership.1   = first(x$ownership),
            teaching.yn.1 = first(x$teaching.yn),
            freestanding.1= first(x$freestanding),
            stringsAsFactors=FALSE
        )
}))
kmain <- kmain[complete.cases(kmain), ]
row.names(kmain) <- NULL

set.seed(1234)
clusters <- kmeans(kmain[2:5], centers=2L)
kmain$cluster <- clusters$cluster

# Swap if the clusters came out backwards
if (sum(kmain$cluster == 1L) == 225L & sum(kmain$cluster == 2L) == 178L) {
    kmain$cluster[kmain$cluster == 1L] <- 0L
    kmain$cluster[kmain$cluster == 2L] <- 1L
    kmain$cluster[kmain$cluster == 0L] <- 2L
}

output$title("Box-and-whisker plots of outcome measures by cluster")

km.plt <- do.call(rbind, lapply(outcome.measures[-c(1,2,5)], function(m) {
    m <- paste0(m, ".avg")
    data.frame(
        measure=m,
        cluster=kmain$cluster,
        value=kmain[[m]],
        stringsAsFactors=FALSE)
}))
km.plt <- km.plt[complete.cases(km.plt), ]


print(lattice::bwplot(
    cluster ~ value | measure,
    data=km.plt,
    auto.key=TRUE,
    scales=list(relation="free"),
    ylab="Cluster",
    xlab="Outcome measure"
))
output$png()


output$title("Means, standard deviations, and two-sample t-tests comparing graded outcome measures between the two clusters identified by the k-means algorithm applied to the graded measures aggregated to the facility")
t.km <- do.call(rbind, lapply(km.outcome.measures, function(u) {
    measure <- kmain[[paste0(u, ".avg")]]
    x <- measure[kmain[["cluster"]] == 1L]
    y <- measure[kmain[["cluster"]] == 2L]
    t <- t.test(y, x)

    result <- array("", dim=c(1L, 4L))
    result[, 1] <- sprintf(
        "%.2f (%.2f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
    result[, 2] <- sprintf(
        "%.2f (%.2f)", mean(y, na.rm=TRUE), sd(y, na.rm=TRUE))
    result[, 3] <- sprintf("%.2f", t$statistic)
    result[, 4] <- sprintf("%s", format.p.value(t$p.value))
    rownames(result) <- u
    colnames(result) <- c("Cluster 1", "Cluster 2", "t-statistic", "p-value")
    result
}))

colnames(t.km)[1] <- sprintf("Cluster 1 (n=%d)", sum(kmain$cluster == 1L))
colnames(t.km)[2] <- sprintf("Cluster 1 (n=%d)", sum(kmain$cluster == 2L))

output$print(t.km)


output$h2("Facility characteristics as predictors of grades")
output$title("Odds ratios of having a <q>Better than national rating or benchmark</q> grade by facility characteristic.")


# Sort main for GEE if we decide to go down that path
gee.main <- main[order(main$facility.id, main$year), ]
gee.main$n.discharge <- gee.main$n.discharge/100
gee.main$population.severity <- gee.main$population.severity*10

t2 <- function(outcome, predictor) {
    fml <- as.formula(paste0(outcome, " ~ ", predictor))

    fit <- try({
        gee::gee(
            fml,
            id=facility.id,
            data=gee.main,
            family=binomial(link="logit"),
            corstr="exchangeable")
    }, silent=TRUE)

    if (inherits(fit, "try-error")) {
        col <- c("Estimate", "Naive S.E.", "Naive z", "Robust S.E.", "Robust z")
        b <- array("", dim=c(1L, 3L))
        rownames(b) <- predictor
        colnames(b) <- c("Odds ratio", "95% CI", "p-value")
        return(b)
    }

    estimate <- coef(summary(fit))[-1, , drop=FALSE]
    b <- estimate[, "Estimate"]
    se <- estimate[, "Robust S.E."]
    z <- estimate[, "Robust z"]

    odds <- exp(b)
    lcl <- exp(b - qnorm(1 - 0.05/2)*se)
    ucl <- exp(b + qnorm(1 - 0.05/2)*se)
    p <- 2*pnorm(-abs(z))

    b <- array("", dim=c(length(b), 3L))
    b[, 1] <- sprintf("%.2f", odds)
    b[, 2] <- sprintf("(%.2f, %.2f)", lcl, ucl)
    b[, 3] <- sprintf("%s", format.p.value(p))

    colnames(b) <- c("Odds ratio", "95% CI", "p-value")
    rownames(b) <- rownames(estimate)
    b
}


# Far too few events to do anything with CAUTI grades
#do.call(rbind, lapply(facility.chrs, function(u) t2("cauti.good", u)))
t2.readmit30 <- do.call(
    rbind,
    lapply(facility.chrs, function(u) t2("readmit30.good", u)))

t2.cdi <- do.call(
    rbind,
    lapply(facility.chrs, function(u) t2("cdi.good", u)))

t2.home <- do.call(
    rbind,
    lapply(facility.chrs, function(u) t2("home.good", u)))

col <- colnames(t2.home)
t2.table <- cbind(t2.readmit30, t2.cdi, t2.home)
colnames(t2.table) <- c(
    paste0("readmit30:", col),
    paste0("cdi:", col),
    paste0("home:", col))

output$print(t2.table)

output$put("<p>Notes:<br/><ol><li>All models are logistic regressions of the outcome with the characteristic as predictor, with GEE being used to perform the estimation with exchangeable working correlation structure. The robust SEs were used for the confidence intervals and p-values.</li><li>CAUTI not included in the table because it did not have enough better than national grades.</li><li>n.discharge is in units of 100, population.severity in units of 0.1.</li><li>Rural not reported for CDI because there were no rural facilities that had a better than national grade for CDI.</li></ol></p>")


#----------------------------------------------------------------------#
# Aim 3: Change over time
#----------------------------------------------------------------------#

output$h1("Change-over-time in relation to facility characteristics")

output$title("Type III ANOVA significance tests of (1) facility characteristics and year in linear mixed effects models of pressure ulcers (v1), CAUTI, readmission within 30 days and CDI, and (2) facility characteristics in linear models of pressure ulcers (v2), rate of falls with injury, and successful return home.")


t4.fit <- function(outcome, predictor, method="lme", type="anova") {

    i <- complete.cases(gee.main[c(outcome, predictor, "year")])
    id <- gee.main[i, "facility.id"]
    y <- gee.main[i, outcome]
    x <- gee.main[i, predictor]

    # Adjust contrasts to get Type III based on
    # https://stats.idre.ucla.edu/r/faq/how-can-i-get-type-iii-tests-of-fixed-effects-in-r
    if (is.factor(x) & type == "anova") contrasts(x) <- contr.sum

    assign(predictor, x)
    time <- gee.main[i, "year"]

    t0 <- min(time)
    Year <- time - t0  # Hope for numerical stability
    degenerate <- (length(unique(Year)) < 2L)


#    if (method == "jrfit") {
#        X <- model.matrix(as.formula(paste0("~ Year*", predictor)))
#        model <- jrfit(X, y, block=id, var.type="sandwich")
#        fit <- summary(model)$coefficients
#    } else if (method == "lme") {
        if (!degenerate) {
            fml <- as.formula(paste0("y ~ Year*", predictor))
            model <- lme(fml, random=~1|id, na.action=na.omit, method="ML")

            # lme needs to be coddled for anova to work
            modelx <- lme(y ~ Year*x, random=~1|id, na.action=na.omit, method="ML", data=data.frame(x=x, Year=Year, y=y, id=id))
            model0 <- lme(y ~ Year, random=~1|id, na.action=na.omit, method="ML")
            modelt <- lme(y ~ x, random=~1|id, na.action=na.omit, method="ML")
            p.value <- anova(model0, modelx)[["p-value"]][2]
            p.trend <- anova(modelt, modelx)[["p-value"]][2]
            col <- c("Estimate", "Std. Error", "DF", "t-value", "p.value")
        } else {
            fml <- as.formula(paste0("y ~ ", predictor))
            model <- lm(fml)
            modelx <- lm(y ~ x, data=data.frame(y=y, x=x))  # To be consistent
            model0 <- lm(y ~ 1)
            p.value <- anova(model0, model)[["Pr(>F)"]][2]
            p.trend <- NULL
            col <- c("Estimate", "Std. Error", "t-value", "p.value")
        }
        fit <- coef(summary(model))
        colnames(fit) <- col
#    } else if (method == "gee") {
#        fml <- as.formula(paste0("y ~ Year*", predictor))
#        model <- gee::gee(fml, id=id, corstr="exchangeable")
#        fit <- coef(summary(model))[-1, , drop=FALSE]
#        fit <- cbind(fit, 2*pnorm(-abs(fit[, "Robust z"])))
#        colnames(fit) <- c("Estimate", "Naive SE", "Naive z", "Std. Error", "Robust z", "p.value")
#    }

    if (type == "coefficients") {
        b <- fit[, "Estimate"]
        se <- fit[, "Std. Error"]
        lcl <- b - qnorm(1 - 0.05/2)*se
        ucl <- b + qnorm(1 - 0.05/2)*se
        p.value <- fit[, "p.value"]

        nsf <- "%.2f"
        s <- rbind(
            c("", "", ""),
            cbind(
                sprintf(nsf, b),
                sprintf("(%.2f, %.2f)", lcl, ucl),
                sprintf("%s", format.p.value(p.value))))
        rownames(s) <- c(predictor, rownames(fit))
        colnames(s) <- paste0(outcome, ":", c("Estimate", "95% CI", "p-value"))
    } else if (type == "anova") {
        if (!degenerate) {
            ftable <- anova(model, type="marginal")
        } else {

            f.aov <- anova(model)
            ftable <- array(0, dim=c(4L, 4L))

            colnames(ftable) <- c("numDF", "denDF", "F-value", "p-value")
            ftable[3L, "numDF"] <- f.aov[["Df"]][1]
            ftable[3L, "denDF"] <- f.aov[["Df"]][2]
            ftable[3L, "F-value"] <- f.aov[["F value"]][1]
            ftable[3L, "p-value"] <- f.aov[["Pr(>F)"]][1]

        }

        s <- array("", dim=c(nrow(ftable) + 1L, 1L))
        p <- format.p.value(ftable[, "p-value"])
        p.compare <- ifelse(grepl("<", p), "", "= ")
        s[-1, 1] <- sprintf(
            "F<sub>%d,%d</sub> = %.2f<br/>p %s%s",
            ftable[, "numDF"], ftable[, "denDF"], ftable[, "F-value"],
            p.compare, p)

        if (degenerate) s[-4L, ] <- ""

        rownames(s) <- c(
            predictor,
            "(Intercept)",
            "Year", "Main effect", "Year x Main effect")
        colnames(s) <- outcome
    }


    if (is.factor(x)) {
        xn <- unique(x)
        if (type == "anova") contrasts(xn) <- contr.sum
    } else {
        xn <- unname(quantile(x, na.rm=TRUE))[c(2L, 4L)]
    }
    p <- expand.grid(x=xn, Year=seq_along(unique(Year))-1L)
    predictions <- augur(modelx, p)
    p$y <- predictions$fit[, 1]
    p$se <- predictions$se.fit
    p$Year <- p$Year + t0

    names(p) <- c(predictor, "year", outcome, "se.fit")


    list(
        summary=s[c(predictor, "(Intercept)", "Main effect", "Year", "Year x Main effect"), , drop=FALSE],
        p.value=p.value,
        p.trend=p.trend,
        model=model,
        prediction=p,
        contrasts=augur.compare(predictions)
    )
}

t4.fits <- lapply(
    as.set(outcome.measures, process.measures),
    function(y) lapply(as.set(facility.chrs), function(x) t4.fit(y, x)))

t4.outcomes <- do.call(cbind,
    lapply(
        outcome.measures,
        function(y) do.call(rbind, lapply(
            t4.fits[[y]],
            function(f) f$summary))
    )
)

output$print(t4.outcomes[rownames(t4.outcomes) != "(Intercept)", , drop=FALSE])


output$title("Type III ANOVA significance tests of (1) facility characteristics and year in linear mixed effects models of patient flu vaccination and health care provider flu vaccination, and (2) facility characteristics in linear models of functional assessment and MSPB.")


t4.process <- do.call(cbind,
    lapply(
        process.measures,
        function(y) do.call(rbind, lapply(
            t4.fits[[y]],
            function(f) f$summary))
    )
)

output$print(t4.process[!(rownames(t4.process) %in% c("(Intercept)", "Year", "Year x Main effect")), , drop=FALSE])


output$h1("Simple effects")

t5.list <- list()

for (measure in outcome.measures) {
    for (chr in facility.chrs) {
        plt.frame <- prediction <- t4.fits[[measure]][[chr]]$prediction
        p.value <- format.p.value(t4.fits[[measure]][[chr]]$p.value)


        if (length(unique(plt.frame$year)) < 2L) next

        if (!(measure %in% names(t5.list))) t5.list[[measure]] <- list()

        colnames(plt.frame) <- c("x", "t", "y", "se")

        plt <-lattice::xyplot(
            y ~ t, group=x, data=plt.frame, type="l",
            auto.key=TRUE,
            scales=list(x=list(at=2016L:2018L)),
            ylab=measure, xlab="Year")

        text <- sprintf("<p>Expected value of %s as a function of year by %s (p %s%s)</p>", measure, chr, ifelse(grepl("<", p.value), "", "= "), gsub("<", "&lt;", p.value))

        t5.list[[measure]][[chr]] <- list(
            means=prediction,
            p.values=t4.fits[[measure]][[chr]]$contrasts,
            plot=plt,
            text=text)

    }
}


# Table of means, SEs and significance results
t5.numbers <- do.call(cbind, lapply(names(t5.list), function(x) {
    do.call(rbind, lapply(names(t5.list[[x]]), function(y) {
        measure.name <- x
        chr.name <- y


        predictions <- t5.list[[x]][[y]]$means
        p.values <- t5.list[[x]][[y]]$p.values$b.values

# Compact letter displays

sig.code <- try({
cld$insert_absorb(setNames(p.values[upper.tri(p.values)] < 0.05, paste0(row(p.values), "-", col(p.values))[upper.tri(p.values)]), lvl_order=order(predictions[, 3]))$Letters
}, silent=TRUE)
if (inherits(sig.code, "try-error")) sig.code <- rep("!", nrow(p.values))

        if (min(predictions[, 2]) == 2016L) {
            row2018 <- predictions[predictions[, 2] == 2016L, , drop=FALSE]
            row2018[, 2] <- 2018L
            row2018[, 3] <- NA
            row2018[, 4] <- NA
            predictions <- rbind(predictions, row2018)
            sig.code <- c(sig.code, rep("", nrow(row2018)))
            for (i in 1:nrow(row2018)) p.values <- rbind(cbind(p.values, 2), 2)
        } else {
            row2016 <- predictions[predictions[, 2] == 2018L, ]
            row2016[, 2] <- 2016L
            row2016[, 3] <- NA
            row2016[, 4] <- NA
            predictions <- rbind(row2016, predictions)
            sig.code <- c(rep("", nrow(row2016)), sig.code)
            for (i in 1:nrow(row2016)) p.values <- rbind(2, cbind(2, p.values))
        }


        i <- do.call(order, predictions[1:2])
        predictions <- predictions[i, ]
        p.values <- p.values[i, i]
sig.code <- sig.code[i]

        ltr <- c(letters, paste0("a", letters))
#        sig.code <- sapply(1:ncol(p.values), function(u) {
#            paste(ltr[which(p.values[, u] < 0.05)], collapse=",")
#        })

        value <- predictions[, 1]
        if (is.factor(value)) {
            value <- as.character(levels(value))[value]
        }

        chr <- as.character(value)
        year <- as.character(predictions[, 2])
        mean <- sprintf("%.2f (%.2f)<sup>%s</sup>", predictions[, 3], predictions[, 4], sig.code)

        result <- cbind(measure.name, chr.name, value, year, ltr[seq_along(mean)], mean)
        colnames(result) <- paste0(x, ".", c("Measure", "Characteristic", "Value", "Year", "Letter", "Mean"))
        result
    }))
}))

t5.columns <- c(
    "bedsores.v1.Characteristic",
    "bedsores.v1.Value",
    "bedsores.v1.Year",
    "bedsores.v1.Letter",
    paste0(c("bedsores.v1", "cauti", "readmit30", "cdi"), ".", "Mean"))

t5 <- t5.numbers[, t5.columns]
t5[t5 == "NA (NA)<sup></sup>"] <- ""

colnames(t5) <- c(
    "Facility.characteristic",
    "Value",
    "Year",
    "Letter",
    "bedsores.v1", "cauti", "readmit30", "cdi")

output$title("Comparison of means from mixed effects models of outcome measures. First and third quartiles are used in the comparison of continuous facility characteristics. P-values are computed using Wald tests of the relevant contrasts. No correction has been made for multiple comparisons.")

output$print(t5)


output$h1("Profile plots")

value.labels <- with(
    read.csv("../data/irfmain_labels.csv", stringsAsFactors=FALSE), {
        setNames(Label, Variable)
})

profile.plot <- function(x) {

    colors <- c("#4472C4", "#ED7D31", "#A5A5A5",
                "#FFC000", "#5B9BD5", "#70AD47",
                "#264478", "#9E480E", "#636363")
    cex <- 1.5
    lwd <- 3.0

    group <- x[, 1]
    if (is.numeric(group)) group <- factor(group, levels=sort(unique(group)), labels=c("Low", "High"))
    year <- x[, 2]
    value <- x[, 3]

    group.name <- value.labels[colnames(x)[1]]
    value.name <- value.labels[colnames(x)[3]]


    plot(year, value, type="n", axes=FALSE, xlab="Year", ylab=value.name, cex.lab=cex)
    axis(1, at=sort(unique(year)), cex.axis=cex)
    axis(2, at=pretty(value, eps.correct=2), cex.axis=cex)

    gi <- 0L
    levels <- unique(group)

    for (g in levels) {
        gi <- gi + 1L
        i <- which(group == g)
        lines(year[i], value[i], lty=gi, col=colors[gi], lwd=lwd)
    }

    legend("bottom", legend=levels,
        title=group.name,
        lty=seq_along(levels), col=colors[seq_along(levels)],
        lwd=lwd, cex=cex, inset=c(0, 1), xpd=TRUE, horiz=TRUE, bty="n")
}


for (measure in names(t5.list)) {
    output$h2(measure)
    for (chr in names(t5.list[[measure]])) {
        output$h3(chr)
        output$put(t5.list[[measure]][[chr]]$text)
        profile.plot(t5.list[[measure]][[chr]][["means"]])
        #print(t5.list[[measure]][[chr]]$plot)
        output$png()
    }
}

for (measure in names(t5.list)) {
    for (chr in names(t5.list[[measure]])) {
        png(sprintf("../output/img/%s_%s.png", measure, chr))
        profile.plot(t5.list[[measure]][[chr]][["means"]])
        dev.off()
    }
}

output$close()