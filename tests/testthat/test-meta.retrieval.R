context("Test: meta.retrieval()")

test_that("The meta retrieval fails when wanted",{

  skip_on_cran()
  skip_on_travis()
  dir <- tempfile()
  expect_error(meta.retrieval(kingdom = "gff",
                              group = "Gammaproteobacteria",
                              db = "gff",
                              type = "assembly_stats",
                              path = dir, combine = TRUE,
                              max_species = 2, mute_citation = TRUE))
  expect_error(meta.retrieval(kingdom = "bacteria",
                          group = "Gammaproteobacteria",
                          db = "gff",
                          type = "assembly_stats",
                          path = dir, combine = TRUE,
                          max_species = 2, mute_citation = TRUE))
  expect_error(meta.retrieval(kingdom = "bacteria",
                              group = "Gammaproteobacteria",
                              db = "refseq",
                              type = "gff",
                              path = dir, combine = TRUE,
                              max_species = 2, mute_citation = TRUE))
})

test_that("The meta retrieval gff subset..",{

  skip_on_cran()
  skip_on_travis()
  dir <- tempfile()
  paths <- meta.retrieval(kingdom = "bacteria",
                          group = "Gammaproteobacteria",
                          db = "refseq",
                          type = "gff",
                          path = dir,
                          max_species = 2, mute_citation = TRUE)
  expect_is(paths, "character")
  expect_equal(length(paths),  2)

  paths2 <- meta.retrieval(kingdom = "bacteria",
                           group = "Gammaproteobacteria",
                           db = "refseq",
                           type = "gff",
                           path = dir,
                           max_species = 2, mute_citation = TRUE)
  expect_is(paths2, "character")
  expect_equal(length(paths2),  2)
  expect_equal(paths,  paths2)

  # Clean retrieval works:
  clean_paths <- clean.retrieval(paths)
  expect_is(paths, "character")
  expect_equal(length(paths),  2)
  expect_equal(tools::file_ext(paths),  c("gz", "gz"))
  expect_equal(tools::file_ext(clean_paths),  c("gff", "gff"))
})

test_that("The meta retrieval assembly_stats subset..",{

  skip_on_cran()
  skip_on_travis()
  dir <- tempfile()
  paths <- meta.retrieval(kingdom = "bacteria",
     group = "Gammaproteobacteria",
     db = "refseq",
     type = "assembly_stats",
     path = dir,
     max_species = 2, mute_citation = TRUE)
  expect_is(paths, "character")
  expect_equal(length(paths),  2)

  paths2 <- meta.retrieval(kingdom = "bacteria",
                          group = "Gammaproteobacteria",
                          db = "refseq",
                          type = "assembly_stats",
                          path = dir,
                          max_species = 2, mute_citation = TRUE)
  expect_is(paths2, "character")
  expect_equal(length(paths2),  2)
  expect_equal(paths,  paths2)
})

test_that("The meta retrieval assembly_stats subset combine..",{

  skip_on_cran()
  skip_on_travis()
  dir <- tempfile()
  paths <- meta.retrieval(kingdom = "bacteria",
                          group = "Gammaproteobacteria",
                          db = "refseq",
                          type = "assembly_stats",
                          path = dir, combine = TRUE,
                          max_species = 2, mute_citation = TRUE)
  expect_is(paths, "data.frame")
  expect_equal(nrow(paths),  2)

  paths2 <- meta.retrieval(kingdom = "bacteria",
                           group = "Gammaproteobacteria",
                           db = "refseq",
                           type = "assembly_stats",
                           path = dir, combine = TRUE,
                           max_species = 2, mute_citation = TRUE)
  expect_is(paths2, "data.frame")
  expect_equal(nrow(paths2),  2)
  expect_equal(paths,  paths2)
})
