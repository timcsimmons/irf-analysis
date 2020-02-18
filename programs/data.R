#----------------------------------------------------------------------#
# Program: data.R
# Author: Tim Simmons
# Date: 2019-04-01
# Purpose: Prepare data for analysis
#----------------------------------------------------------------------#

source("functions.R")

irfanalytic <- read.tbl("../data/irfanalytic.csv")


# Data collection dates for the measures
irfdates <- do.call(rbind, lapply(
    grep("start", names(irfanalytic), value=TRUE),
    function(x) {
        measure <- gsub("\\.start", "", x)
        dates <- unique(
            irfanalytic[c(
                "year",
                paste0(measure, ".start"),
                paste0(measure, ".end"))])
        dates <- dates[complete.cases(dates), , drop=FALSE]
        dates$measure <- measure
        names(dates) <- c("Year", "Start", "End", "Measure")
        dates[c(4, 1:3)]
    })
)

write.csv(irfdates, file="../data/irfdates.csv", row.names=FALSE)


# Gather up notes into a separate data frame
irfnotes <- do.call(rbind, lapply(
    grep("note$", names(irfanalytic), value=TRUE),
    function(x) {
        measure <- gsub("\\.note$", "", x)
        i <- !is.na(irfanalytic[[x]])
        note <- irfanalytic[i, c("CMS_cert_no", "year", x)]
        note$measure <- measure
        names(note) <- c("FaciltyId", "Year", "Note", "Measure")
        note[c(1, 2, 4, 3)]
    })
)

write.csv(irfnotes, file="../data/irfnotes.csv", row.names=FALSE)


# Prepare data frame for main analyses
main <- data.frame(row.names=1:nrow(irfanalytic))

main[["facility.id"]]         <- irfanalytic[["CMS_cert_no"]]
main[["year"]]                <- irfanalytic[["year"]]
main[["facility.name"]]       <- irfanalytic[["Fac_name"]]
main[["city"]]                <- irfanalytic[["City"]]
main[["state"]]               <- irfanalytic[["State"]]
main[["zip"]]                 <- irfanalytic[["Zip"]]
main[["region"]]              <- irfanalytic[["Region"]]
main[["census"]]              <- irfanalytic[["census"]]

main[["n.discharge"]]         <- irfanalytic[["N"]]
main[["cost.charge.ratio"]]   <- irfanalytic[["cost_charge"]]
main[["wage.index"]]          <- irfanalytic[["wage_indx"]]
main[["rural"]]               <- irfanalytic[["Rural"]]
main[["dsh"]]                 <- irfanalytic[["dis_share"]]
main[["ownership"]]           <- irfanalytic[["ownership"]]
main[["teaching"]]            <- irfanalytic[["teaching"]]
main[["teaching.yn"]]         <- as.integer(irfanalytic[["teaching"]] > 0)
main[["freestanding"]]        <- irfanalytic[["freestanding"]]
main[["population.severity"]] <- irfanalytic[["avg_wt"]]


# Patient outcome measures
main[["bedsores.v1"]]         <- irfanalytic[["I_001_01"]]
main[["bedsores.v2"]]         <- irfanalytic[["I_001_02"]]
main[["cauti"]]               <- irfanalytic[["I_006_01"]]
main[["readmit30"]]           <- irfanalytic[["I_007_01"]]
main[["falls"]]               <- irfanalytic[["I_013_01"]]
# (I_014_01 (MRSA) was essentially not reported in any year)
main[["cdi"]]                 <- irfanalytic[["I_015_01"]]
main[["home"]]                <- irfanalytic[["I_019_01"]]


# Process measures
# Assume I_002_01.obs was supposed to be I_002_01
main[["pt.flu"]]              <- irfanalytic[["I_002_01.obs"]]
main[["functional"]]          <- irfanalytic[["I_008_01"]]
main[["hc.flu"]]              <- irfanalytic[["I_016_01"]]
main[["mspb"]]                <- irfanalytic[["I_020_01"]]

main[["cauti.grade"]]         <- irfanalytic[["I_006_01.grade"]]
main[["readmit30.grade"]]     <- irfanalytic[["I_007_01.grade"]]
main[["cdi.grade"]]           <- irfanalytic[["I_015_01.grade"]]
main[["home.grade"]]          <- irfanalytic[["I_019_01.grade"]]

main[["cauti.good"]]          <- is.good.grade(main$cauti.grade)
main[["readmit30.good"]]      <- is.good.grade(main$readmit30.grade)
main[["cdi.good"]]            <- is.good.grade(main$cdi.grade)
main[["home.good"]]           <- is.good.grade(main$home.grade)


# Based on the dates, 2017 values were accidentally carried forward into
# 2018 and so need to be zapped
main[main$year == 2018L, "pt.flu"] <- NA
main[main$year == 2018L, "hc.flu"] <- NA

write.csv(main, file="../data/irfmain.csv", row.names=FALSE)
