# identify_isotopologues() rebuilds the parent from enviPat's NEUTRAL elemental
# mass, while Stage 1 matched the peak using the compound table's
# electron-corrected ion mass (docs/xMSannotator_Input_Formats.md). For ordinary
# adducts the two agree, because adduct_table$mass already carries the electron
# term -- the proton shift is 1.007276, a hydrogen atom less one electron. A
# bare-ion pseudo-adduct has mass = 0, so nothing re-applied the correction and
# the labeller anchored |z| electron masses high: 5.90 ppm for a dication.

ELECTRON <- 0.00054857990924

test_that("ordinary adducts get no electron shift", {
  for (a in c("M+H", "M-H", "M+Na", "M+NH4", "M-H2O-H", "2M+H", "M")) {
    expect_equal(bare_ion_electron_shift(a, 1), 0)
  }
})

test_that("a bare cation is lighter than its neutral formula by |z| electrons", {
  expect_equal(bare_ion_electron_shift("[M]+", 1), -ELECTRON)
  expect_equal(bare_ion_electron_shift("[M]+2", 2), -2 * ELECTRON)
  expect_equal(bare_ion_electron_shift("[M]+3", 3), -3 * ELECTRON)
})

test_that("a bare anion is heavier by the same amount", {
  expect_equal(bare_ion_electron_shift("[M]-", 1), ELECTRON)
  expect_equal(bare_ion_electron_shift("[M]-2", 2), 2 * ELECTRON)
  expect_equal(bare_ion_electron_shift("[M]-5", 5), 5 * ELECTRON)
})

test_that("the shift is taken from the adduct name, not the sign of charge", {
  # adduct_table stores charge unsigned (add_bare_ion_adducts uses abs()), so a
  # positive `charge` must still yield a positive shift for an anion
  expect_equal(bare_ion_electron_shift("[M]-2", 2), 2 * ELECTRON)
  expect_equal(bare_ion_electron_shift("[M]-2", -2), 2 * ELECTRON)
})

# --- the regression: paraquat ------------------------------------------------
# C12H14N2(2+), neutral formula mass 186.11570, ion mass 186.11460.
# Its 13C isotopologue sits at (186.11570 + 1.0033548 - 2e) / 2 = 93.55898.

paraquat_frame <- function() {
  neutral <- 12 * 12 + 14 * 1.00782503 + 2 * 14.0030740
  iso_neutral <- neutral + 1.0033548
  data.frame(
    compound_id = c("INFDPOAKFNIJBF", "INFDPOAKFNIJBF"),
    Formula = c("C12H14N2", "C12H14N2"),
    Adduct = c("[M]+2", "[M]+2_[+1]"),
    mz = c((neutral - 2 * ELECTRON) / 2, (iso_neutral - 2 * ELECTRON) / 2),
    Module_RTclust = c("7_1", "7_1"),
    mean_int_vec = c(1e6, 1.34e5),   # ~13.4%, the 13C abundance for 12 carbons
    stringsAsFactors = FALSE
  )
}

paraquat_adducts <- function() {
  data.frame(adduct = "[M]+2", charge = 2L, factor = 1L, mass = 0,
             stringsAsFactors = FALSE)
}

test_that("a paraquat 13C isotope is labelled at 5 ppm", {
  out <- identify_isotopologues(paraquat_frame(), paraquat_adducts(),
                                isotope_mass_tolerance_ppm = 5)
  expect_false(is.na(out$isotopologue[2]))
  expect_match(out$isotopologue[2], "13C")
})

test_that("the anchor error the fix removes is about 5.9 ppm", {
  # what the labeller computed before: neutral mass, no electron term
  neutral <- 12 * 12 + 14 * 1.00782503 + 2 * 14.0030740
  uncorrected <- (neutral + 1.0033548) / 2
  corrected <- (neutral + 1.0033548 - 2 * ELECTRON) / 2
  expect_equal(1e6 * (uncorrected - corrected) / corrected, 5.9, tolerance = 0.02)
  # and that error alone exceeds a 5 ppm window, so the isotope was rejected
  expect_gt(1e6 * (uncorrected - corrected) / corrected, 5)
})

test_that("neutral compounds are unaffected", {
  # caffeine M+H and its 13C isotope: adduct mass carries the electron term
  d <- data.frame(
    compound_id = c("RYYVLZVUVIJVGH", "RYYVLZVUVIJVGH"),
    Formula = c("C8H10N4O2", "C8H10N4O2"),
    Adduct = c("M+H", "M+H_[+1]"),
    mz = c(195.087652, 196.091007),
    Module_RTclust = c("3_1", "3_1"),
    mean_int_vec = c(1e6, 9.0e4),
    stringsAsFactors = FALSE
  )
  at <- data.frame(adduct = "M+H", charge = 1L, factor = 1L, mass = 1.007276,
                   stringsAsFactors = FALSE)
  out <- identify_isotopologues(d, at, isotope_mass_tolerance_ppm = 5)
  expect_false(is.na(out$isotopologue[2]))
  expect_match(out$isotopologue[2], "13C")
})
