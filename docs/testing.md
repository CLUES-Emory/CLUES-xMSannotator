# Testing

Tests are implemented using the `testthat` and `patrick` framework, which enables parameterized testcases.
Code coverage reports for the various branches of the project are available [here](https://app.codecov.io/gh/CLUES-Emory/CLUES-xMSannotator) and are created with the `covr` package.

Larger testdata is stored on our [GitLab](https://gitlab.ics.muni.cz/umsa/umsa-files/-/tree/master/testdata/recetox-xMSannotator) and downloading of large files required for tests can be automated.

Tests can be run using `devtools::test()` with the working directory in the R project folder. Code coverage can be computed using `covr::package_coverage()` for an inline report. Reports can be uploaded to codecov using `covr::codecov(token = "9fd66806-bd9a-4485-b656-0a0002a38d7e")`.