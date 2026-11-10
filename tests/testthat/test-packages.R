test_that("positionR_packages() lists the four packages in suite-list order", {
  pkgs <- positionR_packages()
  expect_s3_class(pkgs, "positionR_packages")
  expect_s3_class(pkgs, "data.frame")
  expect_named(pkgs, c("package", "role", "version", "site"))
  expect_identical(pkgs$package,
                   c("detqc", "detfield", "arraydesign", "arraydesignR"))
  expect_identical(pkgs$site,
                   paste0("https://jakebrownscombe.github.io/", pkgs$package,
                          "/"))
})

test_that("the roles are the suite list wording", {
  roles <- setNames(positionR_packages()$role, positionR_packages()$package)
  expect_identical(roles[["detqc"]], "read, check and filter detections")
  expect_identical(roles[["detfield"]], "detection fields and positions")
  expect_identical(roles[["arraydesign"]],
                   "design and test the performance of tracking systems")
  expect_identical(roles[["arraydesignR"]],
                   "arraydesign as a point-and-click application")
})

test_that("versions are read from the library and NA when absent", {
  pkgs <- positionR_packages()
  expect_type(pkgs$version, "character")
  # The four are Imports, so the check installs them
  for (p in pkgs$package) {
    expect_identical(pkgs$version[pkgs$package == p],
                     as.character(utils::packageVersion(p)))
  }
  empty_lib <- empty_library()
  absent <- positionR_packages(lib.loc = empty_lib)
  expect_true(all(is.na(absent$version)))
})

test_that("installed_version() gives NA for a package that does not exist", {
  expect_identical(positionR:::installed_version("notapackage.positionR"),
                   NA_character_)
})

test_that("the print names every package, its role and its site", {
  pkgs <- positionR_packages()
  out <- capture.output(res <- print(pkgs))
  expect_identical(res, pkgs)
  expect_match(out[1], "^positionR: 4 packages, 4 installed$")
  for (i in seq_len(nrow(pkgs))) {
    expect_true(any(out == paste(pkgs$package[i], pkgs$version[i])))
    expect_true(any(out == paste0("  ", pkgs$role[i])))
    expect_true(any(out == paste0("  ", pkgs$site[i])))
  }
  expect_true(all(nchar(out) <= 76))
})

test_that("the print says when a package is not installed", {
  absent <- positionR_packages(lib.loc = empty_library())
  out <- capture.output(print(absent))
  expect_match(out[1], "4 packages, 0 installed")
  expect_true(any(out == "detqc not installed"))
})

test_that("loading positionR attaches none of the four packages", {
  expect_false(any(paste0("package:", positionR_packages()$package) %in%
                     search()))
})
