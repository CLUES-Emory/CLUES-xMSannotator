compute_confidence_levels <- function(c,
                                      chemids,
                                      chemscoremat,
                                      filter.by,
                                      max.rt.diff,
                                      adduct_weights,
                                      max_isp,
                                      min_ions_perchem,
                                      adduct_table,
                                      multimer_abundance_check) {
    cur_chemid <- chemids[c]

    curdata <- chemscoremat[which(chemscoremat$compound_id == cur_chemid), ]
    curdata <- curdata[order(curdata$Adduct), ]
    
    bool_check <- 1

    if (!is.na(filter.by[1])) {
        check_adduct <- which(curdata$Adduct %in% filter.by)
        if (length(check_adduct) <= 0) {
            bool_check <- 0
        }
    }

    Confidence <- 0
    if (bool_check == 1) {
        curdata <- get_confidence_stage4(
                        curdata,
                        max.rt.diff,
                        adduct_weights = adduct_weights,
                        filter.by = filter.by,
                        max_isp = max_isp,
                        min_ions_perchem = min_ions_perchem,
                        adduct_table = adduct_table,
                        multimer_abundance_check = multimer_abundance_check
                    )
        if (!is.na(curdata[1, 1])) {
            Confidence <- as.numeric(as.character(curdata[, 1]))
            if (Confidence[1] < 2) {
                if (length(which(curdata$Adduct %in% adduct_weights[which(as.numeric(adduct_weights[, 2]) > 0), 1])) > 0) {
                    if (curdata$score[1] > 10) {
                        mnum <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% curdata$Adduct), 2])))[1]
                        curdata <- curdata[which(curdata$Adduct %in% adduct_weights[which(as.numeric(as.character(adduct_weights[, 2])) >= mnum), 1]), ]
                        Confidence <- 2
                    }
                }
            }
        }
    } else {
        if (length(which(curdata$Adduct %in% adduct_weights[, 1])) > 0) {
            if (curdata$score[1] >= 10) {
                mnum <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% curdata$Adduct), 2])))[1]
                if (length(which(curdata$Adduct %in% filter.by)) > 0) {
                    curdata <- curdata[which(curdata$Adduct %in% filter.by), ]
                    Confidence <- 2
                }
            }
        }
    }

    if (nrow(curdata) > 1) {
        if (curdata$score[1] < 10) {
            if (length(unique(curdata$Adduct)) < 2) {
                Confidence <- 0
            }
        }
    }

    curdata <- cbind(Confidence, curdata)
    curdata <- as.data.frame(curdata)
    curdata <- curdata[, c("Confidence", "compound_id")]
    curdata <- unique(curdata)
    return(curdata)
}

compute_delta_ppm <- function(chemscoremat_with_confidence) {
    # this is fishy but necessary 
    chemscoremat_with_confidence$mz <- as.numeric(as.character(chemscoremat_with_confidence$mz))
    chemscoremat_with_confidence$theoretical.mz <- as.numeric(as.character(chemscoremat_with_confidence$theoretical.mz))
    
    chemscoremat_with_confidence_temp <- chemscoremat_with_confidence[, c("mz", "theoretical.mz")]
    chemscoremat_with_confidence_temp <- apply(chemscoremat_with_confidence_temp, 1, as.numeric)
    chemscoremat_with_confidence_temp <- t(chemscoremat_with_confidence_temp)
    chemscoremat_with_confidence_temp <- as.data.frame(chemscoremat_with_confidence_temp)
    
    delta_ppm <- apply(chemscoremat_with_confidence_temp, 1, function(x) {
        return(10^6 * abs(x[2] - x[1]) / (x[2]))
    })
    delta_ppm <- round(delta_ppm, 2)
    
    chemscoremat_with_confidence <- cbind(chemscoremat_with_confidence[, 1:8], delta_ppm, chemscoremat_with_confidence[, 9:dim(chemscoremat_with_confidence)[2]])
    chemscoremat_with_confidence <- chemscoremat_with_confidence[order(chemscoremat_with_confidence$Confidence, decreasing = TRUE), ]
    return(chemscoremat_with_confidence)
}

boost_confidence_of_IDs <- function(chemscoremat_with_confidence, boostIDs,
                                    max.mz.diff, max.rt.diff, outloc) {
    cnames_boost <- colnames(boostIDs)
    has_mz <- "mz" %in% cnames_boost
    has_time <- "time" %in% cnames_boost

    if (has_mz || has_time) {
        # Proximity-based matching
        good_ind <- sapply(seq_len(nrow(chemscoremat_with_confidence)), function(i) {
            ann_row <- chemscoremat_with_confidence[i, ]

            # Check ID match first
            id_match <- boostIDs$ID == ann_row$compound_id
            if (!any(id_match)) return(FALSE)

            # Filter by mz if present (using fractional/relative tolerance)
            if (has_mz) {
                mz_match <- abs(boostIDs$mz - ann_row$mz) <= max.mz.diff * pmax(abs(boostIDs$mz), abs(ann_row$mz))
                id_match <- id_match & mz_match
            }

            # Filter by time if present (using absolute tolerance in seconds)
            if (has_time) {
                time_match <- abs(boostIDs$time - ann_row$time) <= max.rt.diff
                id_match <- id_match & time_match
            }

            return(any(id_match))
        })

        good_idx <- which(good_ind)
    } else {
        # ID-only matching
        good_idx <- which(chemscoremat_with_confidence$compound_id %in% boostIDs$ID)
    }

    if (length(good_idx) > 0) {
        chemscoremat_with_confidence$Confidence[good_idx] <- 4
        chemscoremat_with_confidence$score[good_idx] <- chemscoremat_with_confidence$score[good_idx] * 100
    }

    return(chemscoremat_with_confidence)
}

#' @importFrom foreach foreach %do% %dopar%
#' @importFrom dplyr left_join
multilevelannotationstep4 <- function(outloc,
                                      chemscoremat,
                                      max.mz.diff = 5,
                                      max.rt.diff = 30,
                                      adduct_weights = NA,
                                      filter.by = NA,
                                      min_ions_perchem = 1,
                                      boostIDs = NA,
                                      boost.mz.diff = NULL,
                                      boost.rt.diff = NULL,
                                      max_isp = 5,
                                      dbAllinf = NA,
                                      mz_rt_feature_id_map = NULL,
                                      adduct_table = NULL,
                                      multimer_abundance_check = TRUE) {
    chemids <- unique(chemscoremat$compound_id)

    # Load package adduct_table if not provided
    if (is.null(adduct_table)) {
        data("adduct_table", envir = environment())
    }

    adduct_weights <- create_adduct_weights(adduct_weights)

    # assign confidence level
    chemscoremat_conf_levels <- lapply(
        seq_len(length(chemids)),
        compute_confidence_levels,
        chemids,
        chemscoremat,
        filter.by,
        max.rt.diff,
        adduct_weights,
        max_isp,
        min_ions_perchem,
        adduct_table,
        multimer_abundance_check
    )
    
    chemscoremat_conf_levels <- plyr::ldply(chemscoremat_conf_levels, rbind)
    chemscoremat_conf_levels <- as.data.frame(chemscoremat_conf_levels)

    chemscoremat_with_confidence <- merge(chemscoremat_conf_levels, unique(chemscoremat), by = "compound_id")
    chemscoremat_with_confidence <- as.data.frame(chemscoremat_with_confidence)
    
    cnames1 <- colnames(chemscoremat_with_confidence)
    
    chemscoremat_with_confidence <- compute_delta_ppm(chemscoremat_with_confidence)

    # (II) Presence of required adducts/forms specified by the user
    # for assignment to high confidence categories (e.g., M + H).
    if (!identical(boostIDs, NA)) {
        boost_mz <- if (!is.null(boost.mz.diff)) boost.mz.diff else max.mz.diff
        boost_rt <- if (!is.null(boost.rt.diff)) boost.rt.diff else max.rt.diff
        chemscoremat_with_confidence <- boost_confidence_of_IDs(
            chemscoremat_with_confidence, boostIDs, boost_mz, boost_rt, outloc
        )
    }

    t2 <- table(chemscoremat_with_confidence$mz)
    uniquemz <- names(which(t2 == 1))

    # assign match category
    chemscoremat_with_confidence$MatchCategory <- rep("Multiple", dim(chemscoremat_with_confidence)[1])
    chemscoremat_with_confidence$MatchCategory[which(chemscoremat_with_confidence$mz %in% uniquemz)] <- "Unique"

    # Join feature ID column if mapping provided
    if (!is.null(mz_rt_feature_id_map)) {
        chemscoremat_with_confidence <- left_join(chemscoremat_with_confidence, mz_rt_feature_id_map, by = c("mz", "time"))
    }

    write.table(chemscoremat_with_confidence, file = file.path(outloc, "Stage4_confidence_levels.txt"), sep = "\t", row.names = FALSE)

    chemscoremat_with_confidence <- as.data.frame(chemscoremat_with_confidence)
    chemscoremat_with_confidence <- chemscoremat_with_confidence[order(chemscoremat_with_confidence$Confidence, decreasing = TRUE), ]

    print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
    print(table(chemscoremat_with_confidence$Confidence[-which(duplicated(chemscoremat_with_confidence$compound_id) == TRUE)]))

    print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
    print(table(chemscoremat_with_confidence$Confidence[-which(duplicated(chemscoremat_with_confidence$Formula) == TRUE)]))

    return(chemscoremat_with_confidence)
}