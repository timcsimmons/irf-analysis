#----------------------------------------------------------------------#
# Program: profile.R
# Author: Tim Simmons
# Date: 2019-02-14
# Purpose: Profile data for Graham project
#----------------------------------------------------------------------#

source("functions.R")
source("ods.R")

output <- html("../output/profile.html")

irf <- read.main()

output$h1("Missingness report")
output$title("Rate of observed data data by year")
(function() {
    irf.missing <- irf
    for (col in names(irf.missing)) {
        irf.missing[[col]] <-
            factor(as.integer(!is.na(irf.missing[[col]])), levels=0:1)
    }
    irf.missing$year <- irf$year
    t <- tabulator(irf.missing, ~ year)
    t$n()
    for (col in names(irf.missing)) {
        if (col != "year") t$table(col)
    }

    output$print(t)
})()

output$h1("Descriptive statistics")
output$title("Descriptive summaries of facility characteristics and performance measures by year")
t0 <- tabulator(irf, ~ -1 + year)
t0$n()
t0$mean("cost.charge.ratio", mean="%.2f", sd="%.2f")
t0$mean("wage.index", mean="%.2f", sd="%.2f")
t0$table("rural")
t0$median("dsh")
t0$table("census")
t0$table("ownership")
t0$median("teaching", median="%.2f", q1="%.2f", q3="%.2f")
t0$table("freestanding")
t0$mean("population.severity", mean="%.2f", sd="%.2f")

t0$median("bedsores.v1")
t0$median("bedsores.v2")
t0$median("pt.flu")
t0$median("cauti")
t0$median("readmit30")
t0$median("functional")
t0$median("falls")
t0$median("cdi")
t0$median("hc.flu")
t0$median("home")
t0$median("mspb")
t0$table("cauti.grade")
t0$table("readmit30.grade")
t0$table("cdi.grade")
t0$table("home.grade")

output$print(t0)

for (tbl in t0$result()) {
    col <- row.names(tbl)[1]
    if (!is.null(col)) {
        col <- strsplit(row.names(tbl)[1], "/")[[1]][1]
        if (class(tbl) %in% c("tabulator.mean", "tabulator.median")) hist(irf[[col]], main=col)
        if (class(tbl) == "tabulator.table") barplot(table(irf[[col]]), main=col)
        output$png()
    }
}


output$close()