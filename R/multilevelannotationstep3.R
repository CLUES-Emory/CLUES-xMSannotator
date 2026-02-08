p_test <- function(data) {
  counts <- matrix(
    data = data,
    nrow = 2
  )

  result <- fisher.test(counts)
  return(result$p.value)
}


filter_score_and_adducts <- function(chemscoremat, scorethresh, adduct_weights) {
  indices <- which(
    chemscoremat$score >= scorethresh &
      chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
  )
  return(chemscoremat[indices, ])
}


count_chemicals_occurence <- function(module_data, pathway_chemicals, scorethresh, adduct_weights) {
  indices_in_pathway <- which(module_data$compound_id %in% pathway_chemicals)
  module_data_pathway <- module_data[indices_in_pathway, ]

  module_data <- filter_score_and_adducts(module_data, scorethresh, adduct_weights)
  module_data_pathway <- filter_score_and_adducts(module_data_pathway, scorethresh, adduct_weights)

  num_pathway <- length(unique(module_data_pathway$compound_id))
  num <- length(unique(module_data$compound_id)) - num_pathway
  return(list(num_pathway = num_pathway, num = num))
}


compute_score_pathways <- function(chemscoremat, db, pathwaycheckmode, scorethresh, adduct_weights, max_diff_rt, bad_path_IDs, db_name = "HMDB") {
  pthresh <- 0.05

  pathway_ids <- unique(as.character(db[, 2]))

  module_num <- gsub(chemscoremat$Module_RTclust,
    pattern = "_[0-9]*",
    replacement = ""
  )

  chemscoremat <- cbind(chemscoremat, module_num)

  if (!is.na(pathwaycheckmode)) {
    for (path_id in pathway_ids)
    {
      if (!path_id %in% bad_path_IDs) {
        pathway_chemicals <- db[which(db[, 2] %in% path_id), 1]
        # total number of chemicals in pathway
        all_cur_path_numchem <- length(unique(pathway_chemicals))
        curmchemical <- filter_score_and_adducts(chemscoremat, scorethresh, adduct_weights)

        # chemicals in pathway
        curmchemical_in_pathway <- curmchemical[which(curmchemical$compound_id %in% pathway_chemicals), ]

        # molecules of interest in pathway (a)
        num_chems_inpath <- length(unique(curmchemical_in_pathway$compound_id))

        # non-focus molecules associated with pathway (c)
        num_chem_inpath_notinterest <- all_cur_path_numchem - num_chems_inpath

        # chemicals NOT in pathway
        curmchemical_not_in_pathway <- curmchemical[-which(curmchemical$compound_id %in% pathway_chemicals), ]

        # focus molecules not associated with this pathway (b)
        num_chems_notinpath <- length(unique(curmchemical_not_in_pathway$compound_id))

        all_notcurpath_numchem <- length(db[-which(db[, 1] %in% curmchemical_not_in_pathway$compound_id), 1])

        p_value <- p_test(c(
          num_chems_inpath,
          num_chem_inpath_notinterest,
          num_chems_notinpath,
          all_notcurpath_numchem
        ))

        if (p_value <= pthresh) {
          pathway_chemicals_to_iterate <- pathway_chemicals


          for (chemname in pathway_chemicals_to_iterate) {
            pathway_indices <- which(as.character(curmchemical_in_pathway$compound_id) == chemname)
            curmchemicaldata <- curmchemical_in_pathway[pathway_indices, ]

            t2 <- table(curmchemicaldata$module_num)
            cur_module <- names(t2[which(t2 == max(t2)[1])])

            if (nrow(curmchemicaldata) > 0) {
              if (pathwaycheckmode == "pm") {
                t1 <- table(curmchemical_in_pathway$module_num)
                num_chems <- t1[as.character(cur_module)]
              }

              # if pathwaycheckmode is NOT "pm" this will be an error !!!
              num_chems <- round(num_chems, 0)

              # a and b
              module_indices <- which(chemscoremat$module_num == cur_module)
              cur_module_data <- chemscoremat[module_indices, ]
              result <- count_chemicals_occurence(cur_module_data, pathway_chemicals, scorethresh, adduct_weights)

              cur_module_pathway_number <- result$num_pathway
              cur_module_number <- result$num

              # c and d
              other_module_indices <- which(chemscoremat$module_num != cur_module)
              other_module_data <- chemscoremat[other_module_indices, ]
              result <- count_chemicals_occurence(other_module_data, pathway_chemicals, scorethresh, adduct_weights)

              other_module_pathway_number <- result$num_pathway
              other_module_number <- result$num

              if (cur_module_pathway_number > 1) {
                p_value <- p_test(c(
                  cur_module_pathway_number,
                  other_module_pathway_number,
                  cur_module_number,
                  other_module_number
                ))
              } else {
                p_value <- 1
              }

              if (p_value <= 0.2) {
                if (num_chems >= 3) {
                  if (is.na(curmchemicaldata$score[1])) {
                    diff_rt <- max(curmchemicaldata$time) - min(curmchemicaldata$time)

                    if (diff_rt > max_diff_rt & length(which(t2 > 1)) == 1) {
                      curmchemicaldata$score <- rep(0.1, length(curmchemicaldata$score))
                    } else {
                      curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                    }
                  }
                  chemical_name <- chemname                  

                  if (curmchemicaldata$score[1] < scorethresh) {
                    indices <- which(
                      as.character(chemscoremat$compound_id) == chemical_name &
                        chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                    )
                  } else {
                    indices <- which(as.character(chemscoremat$compound_id) == chemical_name)
                  }

                  chemscoremat$score[indices] <- as.numeric(chemscoremat$score[indices][1]) + num_chems
                }
              }
            }
          }
        }
      }
    }
  }

  chemscoremat <- replace_x_names(chemscoremat)
  return(chemscoremat)
}


preprocess_db <- function(db, index, pattern) {
  temp_matrix <- apply(db, 1, function(x) {
    chemid <- x[1]
    match <- gregexpr(x[index], pattern = pattern)
    regexp_check <- attr(match[[1]], "match.length")
    if (regexp_check[1] < 0) {
      pathid <- "-"
    } else {
      pathid <- strsplit(x = x[index], split = ";")
      pathid <- unlist(pathid)
    }
    return(cbind(chemid, pathid))
  })

  matrix <- plyr::ldply(temp_matrix, rbind)
  return(matrix)
}


select_y_names <- function(chemscoremat, column_names) {
  # y because we want chemCompMZ ID and Name
  column_names_y <- column_names
  column_names_y[1] <- "cur_chem_score"
  column_names_y[7:8] <- paste(column_names_y[7:8], ".y", sep = "")
  column_names_y[11] <- "tempadduct"
  chemscoremat <- chemscoremat[, column_names_y]
  return(chemscoremat)
}


replace_x_names <- function(chemscoremat) {
  cnames <- colnames(chemscoremat)
  cnames <- gsub(cnames, pattern = ".x", replacement = "")
  colnames(chemscoremat) <- cnames
  return(chemscoremat)
}


sanitize_chemscoremat <- function(chemscoremat, chemCompMZ, column_names) {
  chemscoremat$Formula <- gsub(chemscoremat$Formula,
    pattern = "_.*",
    replacement = ""
  )

  tempadduct <- chemscoremat$Adduct

  chemscoremat$Adduct <- gsub(chemscoremat$Adduct,
    pattern = "_.*",
    replacement = ""
  )

  chemscoremat <- cbind(chemscoremat, tempadduct)

  chemscoremat <- merge(chemscoremat,
    chemCompMZ[, c(2:4, 6)],
    by = c("Formula", "Adduct")
  )

  chemscoremat <- select_y_names(chemscoremat, column_names)
  colnames(chemscoremat) <- column_names
  return(chemscoremat)
}


#' @importFrom dplyr left_join
multilevelannotationstep3 <- function(chemscoremat,
                                      adduct_weights = NA,
                                      db_name = "HMDB",
                                      max_diff_rt,
                                      pathwaycheckmode = "p",
                                      scorethresh = 0.1,
                                      outloc = tempdir(),
                                      mz_rt_feature_id_map = NULL,
                                      chemCompMZ = NULL,
                                      pathway_data = NULL,
                                      excluded_pathways = NULL,
                                      excluded_pathway_compounds = NULL) {
  adduct_weights <- create_adduct_weights(adduct_weights)

  column_names <- c(
    "score",
    "Module_RTclust",
    "mz",
    "time",
    "MatchCategory",
    "theoretical.mz",
    "compound_id",
    "Name",
    "Formula",
    "MonoisotopicMass",
    "Adduct",
    "mean_int_vec",
    "MD"
  )

  if (db_name == "HMDB") {
    # HMDB mode: existing logic
    chemscoremat <- sanitize_chemscoremat(chemscoremat, chemCompMZ, column_names)

    hmdbbad <- c("HMDB29244", "HMDB29245", "HMDB29246")
    bad_indices <- which(chemscoremat$compound_id %in% hmdbbad)

    if (length(bad_indices) > 0) {
      chemscoremat <- chemscoremat[-bad_indices, ]
    }

    data(hmdbAllinf)
    hmdbAllinfv3.5 <- hmdbAllinf[, -c(26:27)]
    rm(hmdbAllinf, envir = .GlobalEnv)
    db <- preprocess_db(hmdbAllinfv3.5, 14, "SM")
    chemscoremat <- merge(chemscoremat,
      hmdbAllinfv3.5,
      by.x = "compound_id",
      by.y = "HMDBID"
    )
    rm(hmdbAllinfv3.5)
    bad_path_IDs <- c("-")
  } else if (db_name == "custom") {
    # Custom mode: use user-provided pathway_data
    names(chemscoremat)[names(chemscoremat) == "cur_chem_score"] <- "score"

    # Validate and filter pathway_data
    pathway_data <- as_pathway_table(pathway_data)
    if (!is.null(excluded_pathways)) {
      pathway_data <- pathway_data[!pathway_data$pathway %in% excluded_pathways, ]
    }
    if (!is.null(excluded_pathway_compounds)) {
      pathway_data <- pathway_data[!pathway_data$compound %in% excluded_pathway_compounds, ]
    }

    # Convert to db format: 2-column matrix (compound_id, pathway_id)
    db <- as.matrix(pathway_data[, c("compound", "pathway")])
    bad_path_IDs <- character(0)
  } else {
    stop("db_name must be 'HMDB' or 'custom'")
  }

  chemscoremat <- compute_score_pathways(chemscoremat, db, pathwaycheckmode,
                                          scorethresh, adduct_weights, max_diff_rt,
                                          bad_path_IDs, db_name)

  chemscoremat <- replace_x_names(chemscoremat)

  good_ind <- which(chemscoremat$score >= scorethresh)
  if (length(good_ind) == 0) {
    chemscoremat <- {}
  } else {
    # rotate column names
    column_names[1:7] <- c(column_names[7], column_names[1:6])
    chemscoremat <- chemscoremat[, column_names]
  }

  # Join feature ID column if mapping provided
  if (!is.null(mz_rt_feature_id_map) && length(chemscoremat) > 0) {
    chemscoremat <- left_join(chemscoremat, mz_rt_feature_id_map, by = c("mz", "time"))
  }

  output_file <- if (db_name == "custom") "Stage3_custom_pathways.txt" else "Stage3_HMDB_pathways.txt"
  write.table(chemscoremat, file = file.path(outloc, output_file), sep = "\t", row.names = FALSE)

  return(chemscoremat)
}