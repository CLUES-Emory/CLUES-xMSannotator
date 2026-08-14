test_that("bare_ion_adduct names encode sign and magnitude", {
  expect_equal(
    bare_ion_adduct(c(0L, 1L, 2L, -1L, -5L)),
    c(NA_character_, "[M]+", "[M]+2", "[M]-", "[M]-5")
  )
})

test_that("bare ion names avoid the patterns the pipeline matches on", {
  names <- bare_ion_adduct(c(1L, 2L, -1L, -5L))

  # multilevelannotationstep3 strips everything after an underscore
  expect_false(any(grepl("_", names, fixed = TRUE)))
  # multilevelannotationstep4 multimer test
  expect_false(any(grepl("([2-3]+M)", names)))
  # isotope suffix regex must not fire on the bare name, but must strip cleanly
  isotope_pattern <- "(_\\[(\\+|\\-)[0-9]+\\])"
  expect_false(any(grepl(isotope_pattern, names)))
  expect_equal(gsub(isotope_pattern, "", paste0(names, "_[+1]")), names)
  # must not collide with shipped adduct names
  expect_false(any(names %in% get_default_adduct_table()$adduct))
})

test_that("forms_valid_charge_pair restricts charged compounds to their bare ion", {
  # neutral compounds take real adducts but never a zero-mass-shift adduct
  expect_true(forms_valid_charge_pair(0L, "M+H"))
  expect_false(forms_valid_charge_pair(0L, "M"))
  expect_false(forms_valid_charge_pair(0L, "[M]+2"))

  # charged compounds take only their own bare ion
  expect_true(forms_valid_charge_pair(1L, "[M]+"))
  expect_true(forms_valid_charge_pair(2L, "[M]+2"))
  expect_true(forms_valid_charge_pair(-5L, "[M]-5"))
  expect_false(forms_valid_charge_pair(1L, "M+H"))
  expect_false(forms_valid_charge_pair(1L, "M"))
  expect_false(forms_valid_charge_pair(2L, "[M]+"))
})

test_that("forms_valid_charge_pair is vectorized and never returns NA", {
  result <- forms_valid_charge_pair(
    c(0L, 1L, 2L, -5L, 0L),
    c("M+H", "[M]+", "M", "[M]-5", "[M]+")
  )
  expect_equal(result, c(TRUE, TRUE, FALSE, TRUE, FALSE))
  expect_false(anyNA(result))
})

test_that("as_compound_table defaults charge to 0 and coerces to integer", {
  base <- data.frame(
    compound = 1:3,
    name = c("a", "b", "c"),
    molecular_formula = c("C6H6", "C6H6", "C6H6"),
    monoisotopic_mass = c(78, 78, 78),
    stringsAsFactors = FALSE
  )

  expect_equal(as_compound_table(base)$charge, rep(0L, 3))

  doubles <- transform(base, charge = c(0, 2, -5))
  expect_identical(as_compound_table(doubles)$charge, c(0L, 2L, -5L))

  with_na <- transform(base, charge = c(0L, NA_integer_, 2L))
  expect_identical(as_compound_table(with_na)$charge, c(0L, 0L, 2L))

  expect_error(as_compound_table(transform(base, charge = c("a", "b", "c"))))
})

test_that("as_compound_table is idempotent", {
  data <- data.frame(
    compound = 1:2,
    name = c("a", "b"),
    molecular_formula = c("C6H6", "C12H14N2+2"),
    monoisotopic_mass = c(78, 186.1146),
    charge = c(0L, 2L),
    stringsAsFactors = FALSE
  )
  once <- as_compound_table(data)
  expect_equal(as_compound_table(once), once)
})

test_that("add_bare_ion_adducts adds one row per distinct non-zero charge", {
  adducts <- as_adduct_table(get_default_adduct_table())
  neutral_only <- data.frame(charge = c(0L, 0L))
  no_charge_col <- data.frame(compound = 1:2)

  expect_equal(add_bare_ion_adducts(adducts, neutral_only), adducts)
  expect_equal(add_bare_ion_adducts(adducts, no_charge_col), adducts)

  mixed <- data.frame(charge = c(0L, 2L, 2L, -5L))
  extended <- add_bare_ion_adducts(adducts, mixed)
  expect_equal(nrow(extended), nrow(adducts) + 2)
  expect_setequal(setdiff(extended$adduct, adducts$adduct), c("[M]+2", "[M]-5"))

  added <- extended[extended$adduct == "[M]+2", ]
  expect_equal(added$charge, 2L)
  expect_equal(added$factor, 1L)
  expect_equal(added$mass, 0)

  # idempotent, and the synthesized rows pass validation
  expect_equal(add_bare_ion_adducts(extended, mixed), extended)
  expect_silent(as_adduct_table(extended))
})

test_that("add_bare_ion_adducts accepts charge as character or double", {
  adducts <- as_adduct_table(get_default_adduct_table())
  from_text <- add_bare_ion_adducts(adducts, data.frame(charge = c("0", "+2")))
  from_dbl <- add_bare_ion_adducts(adducts, data.frame(charge = c(0, 2)))
  expect_equal(setdiff(from_text$adduct, adducts$adduct), "[M]+2")
  expect_equal(from_text, from_dbl)
})

test_that("charged compounds are matched only at mass / |z|", {
  # Paraquat is a dication: monoisotopic_mass is the ion mass, so the real
  # observation is at 186.11460 / 2. A neutral of the same mass is the control.
  compounds <- data.frame(
    compound = 1:2,
    name = c("paraquat", "neutral twin"),
    molecular_formula = c("C12H14N2+2", "C12H14N2"),
    monoisotopic_mass = c(186.11460, 186.11460),
    charge = c(2L, 0L),
    stringsAsFactors = FALSE
  )
  peaks <- data.frame(
    mz = c(93.05730, 187.12188, 186.11460),
    rt = c(100, 100, 100)
  )
  adducts <- get_default_adduct_table()
  adducts <- adducts[adducts$adduct %in% c("M", "M+H", "M+2H", "M+Na"), ]

  result <- simple_annotation(peaks, compounds, adducts)

  paraquat <- result[result$name == "paraquat", ]
  expect_equal(nrow(paraquat), 1)
  expect_equal(paraquat$adduct, "[M]+2")
  expect_equal(paraquat$mz, 93.05730)

  twin <- result[result$name == "neutral twin", ]
  expect_true("M+H" %in% twin$adduct)
  # never at its own neutral mass, and never under a bare-ion name
  expect_false(186.11460 %in% twin$mz)
  expect_false(any(startsWith(twin$adduct, "[M]")))

  # charge is an input constraint, not an output column
  expect_false("charge" %in% names(result))
})

test_that("a database without a charge column annotates exactly as before", {
  compounds <- data.frame(
    compound = 1L,
    name = "neutral",
    molecular_formula = "C12H14N2",
    monoisotopic_mass = 186.11460,
    stringsAsFactors = FALSE
  )
  peaks <- data.frame(mz = c(187.12188, 186.11460), rt = c(100, 100))
  adducts <- get_default_adduct_table()
  adducts <- adducts[adducts$adduct %in% c("M+H", "M+Na"), ]

  result <- simple_annotation(peaks, compounds, adducts)
  expect_equal(result$adduct, "M+H")
  expect_false("charge" %in% names(result))
})

test_that("isotope spacing is scaled by the charge of the ion", {
  # M, M+1 and M+2 for a dication sit 0.5017 m/z apart, not 1.0034
  base <- 93.05730
  make_peaks <- function(z) {
    data.frame(
      peak = 1:3,
      mz = c(base, base + 1.003355 / z, base + 2.006710 / z),
      rt = rep(100, 3),
      rt_cluster = rep(1, 3),
      mass_defect = rep(1, 3),
      mean_intensity = c(100, 10, 1)
    )
  }
  make_query <- function(z) {
    dplyr::tibble(
      peak = 0L, expected_mass = base, rt = 100, rt_cluster = 1,
      mass_defect = 1, charge = z, compound = 1, compound_id = "X",
      adduct = "A", molecular_formula = "C12H14N2", name = "x",
      monoisotopic_mass = base, mean_intensity = 100
    )
  }

  for (z in c(1, 2)) {
    found <- filter_isotopes(
      make_query(z), 0.1, make_peaks(z),
      mass_defect_tolerance = 0.1, rt_tolerance = 5
    )
    expect_equal(sort(found$mass_number_difference), c(0, 1, 2))
  }

  # Treated as singly charged, the dication's true M+1 peak (+0.5017) is lost:
  # only the +1.0034 peak survives, and it is mislabelled as M+1.
  wrong <- filter_isotopes(
    make_query(1), 0.1, make_peaks(2),
    mass_defect_tolerance = 0.1, rt_tolerance = 5
  )
  expect_false(any(abs(wrong$mz - (base + 1.003355 / 2)) < 1e-6))
  expect_equal(sort(wrong$mass_number_difference), c(0, 1))
})

test_that("compute_isotopes derives ion charge from the adduct table", {
  annotation <- dplyr::tibble(
    peak = 1L, adduct = "M+2H", compound = 1, expected_mass = 100,
    mz = 100, rt = 50, rt_cluster = 1, mass_defect = 1,
    molecular_formula = "C12H14N2", monoisotopic_mass = 198,
    compound_id = "X", name = "x", mean_intensity = 100
  )
  # 12.98 is the natural M+1 abundance of C12H14N2, so the intensity-ratio
  # filter passes and the test isolates the charge behaviour.
  peaks <- dplyr::tibble(
    peak = 2L, mz = 100 + 1.003355 / 2, rt = 50, rt_cluster = 1,
    mass_defect = 1, mean_intensity = 12.98
  )
  weights <- data.frame(adduct = "M+2H", weight = 5)
  adducts <- get_default_adduct_table()

  with_charge <- compute_isotopes(
    annotation, weights, peak_table = peaks, rt_tolerance = 5,
    mass_defect_tolerance = 0.1, isotope_mass_tolerance_ppm = 5,
    adduct_table = adducts
  )
  expect_true(any(with_charge$mass_number_difference == 1))

  # without an adduct table every ion is assumed singly charged, as before
  without <- compute_isotopes(
    annotation, weights, peak_table = peaks, rt_tolerance = 5,
    mass_defect_tolerance = 0.1, isotope_mass_tolerance_ppm = 5
  )
  expect_false(any(without$mass_number_difference == 1))
})

test_that("step 4 does not duplicate the M adduct row", {
  adduct_table <- as_adduct_table(get_default_adduct_table())
  expect_true("M" %in% adduct_table$adduct)

  monoisot <- data.frame(
    adduct = "M", charge = 1, factor = 1, mass = 0, stringsAsFactors = FALSE
  )
  local_table <- if (!is.element("M", adduct_table$adduct)) {
    rbind(adduct_table, monoisot)
  } else {
    adduct_table
  }
  expect_equal(sum(local_table$adduct == "M"), 1)
})

test_that("get_isotopologue_labels handles charge-suffixed formulas", {
  utils::data(isotopes, package = "enviPat", envir = environment())
  expect_error(get_isotopologue_labels("C12H14N2+2", isotopes), NA)
  expect_equal(
    get_isotopologue_labels("C12H14N2+2", isotopes),
    get_isotopologue_labels("C12H14N2", isotopes)
  )
})


# --- polarity filtering -------------------------------------------------------

mixed_charge_table <- function() {
  data.frame(
    compound = 1:5,
    name = c("neutral", "cation", "dication", "anion", "penta-anion"),
    molecular_formula = rep("C6H6", 5),
    monoisotopic_mass = rep(78, 5),
    charge = c(0L, 1L, 2L, -1L, -5L),
    stringsAsFactors = FALSE
  )
}

test_that("filter_compounds_by_polarity keeps neutrals plus the matching sign", {
  db <- mixed_charge_table()

  pos <- filter_compounds_by_polarity(db, "pos")
  expect_equal(pos$charge, c(0L, 1L, 2L))

  neg <- filter_compounds_by_polarity(db, "neg")
  expect_equal(neg$charge, c(0L, -1L, -5L))
})

test_that("filter_compounds_by_polarity is a no-op without a target polarity", {
  db <- mixed_charge_table()
  expect_equal(suppressWarnings(filter_compounds_by_polarity(db, NULL)), db)

  # a table with no charge column cannot be filtered
  no_charge <- db[, setdiff(names(db), "charge")]
  expect_equal(filter_compounds_by_polarity(no_charge, "pos"), no_charge)
})

test_that("filter_compounds_by_polarity is idempotent and coerces charge", {
  db <- mixed_charge_table()
  once <- filter_compounds_by_polarity(db, "pos")
  expect_equal(filter_compounds_by_polarity(once, "pos"), once)

  as_text <- transform(db, charge = c("0", "+1", "+2", "-1", "-5"))
  as_dbl <- transform(db, charge = c(0, 1, 2, -1, -5))
  expect_equal(nrow(filter_compounds_by_polarity(as_text, "pos")), 3)
  expect_equal(nrow(filter_compounds_by_polarity(as_dbl, "pos")), 3)

  # NA charge is treated as neutral, so it survives either polarity
  with_na <- transform(db, charge = c(NA_integer_, 1L, 2L, -1L, -5L))
  expect_equal(nrow(filter_compounds_by_polarity(with_na, "neg")), 3)
})

test_that("a mixed-sign table warns only when no polarity is given", {
  db <- mixed_charge_table()

  expect_warning(filter_compounds_by_polarity(db, NULL), "charge_type")
  expect_silent(filter_compounds_by_polarity(db, "pos"))

  # single-sign and neutral-only tables have nothing to disambiguate
  expect_silent(filter_compounds_by_polarity(db[db$charge >= 0, ], NULL))
  expect_silent(filter_compounds_by_polarity(db[db$charge == 0, ], NULL))
  expect_silent(filter_compounds_by_polarity(db[, names(db) != "charge"], NULL))
})

test_that("advanced_annotation rejects an unknown charge_type", {
  expect_error(
    advanced_annotation(
      peak_table = data.frame(mz = 100, rt = 1, s1 = 1),
      compound_table = mixed_charge_table(),
      charge_type = "positive"
    ),
    "charge_type"
  )
})

test_that("the opposite polarity is absent from a filtered annotation", {
  compounds <- data.frame(
    compound = 1:3,
    name = c("neutral twin", "paraquat", "an anion"),
    molecular_formula = c("C12H14N2", "C12H14N2+2", "C11H16N3O2-"),
    monoisotopic_mass = c(186.11460, 186.11460, 222.12480),
    charge = c(0L, 2L, -1L),
    stringsAsFactors = FALSE
  )
  peaks <- data.frame(
    mz = c(187.12188, 93.05730, 222.12480),
    rt = c(100, 100, 100)
  )
  adducts <- get_default_adduct_table()
  adducts <- adducts[adducts$adduct %in% c("M+H", "M+Na"), ]

  pos_only <- filter_compounds_by_polarity(compounds, "pos")
  adducts_pos <- add_bare_ion_adducts(adducts, pos_only)
  expect_true(all(c("[M]+2") %in% adducts_pos$adduct))
  expect_false(any(startsWith(adducts_pos$adduct, "[M]-")))

  result <- simple_annotation(peaks, pos_only, adducts_pos)
  expect_false("an anion" %in% result$name)
  expect_true("paraquat" %in% result$name)
  expect_true("neutral twin" %in% result$name)
  expect_false(any(startsWith(result$adduct, "[M]-")))
})
