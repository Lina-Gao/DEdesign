library(testthat)
library(DEdesign)
library(dplyr)

context("design")

test_that("treatment4 replicate4 Latin square",{
  des <- gendesign(treatments= data.frame(trt = letters[1:4], replicates = rep(4,4)), seed=1)
  expect_true(all.equal(
    des$Design$design %>% dplyr::filter(lane==1) %>% dplyr::arrange(adapter) %>% .$trt %>% as.character(),
    c("a", "b", "c", "d")))
  expect_true(all.equal(
    des$Design$design %>% dplyr::filter(lane==2) %>% dplyr::arrange(adapter) %>% .$trt %>% as.character(),
    c("d", "c", "b", "a")))
  expect_true(all.equal(
    des$Design$BlocksEfficiency$`D-efficiency` %>% as.character %>% as.numeric,
    c(1,1)))
})

test_that("treatment6 replicate3",{
  des <- gendesign(treatments=data.frame(trt = letters[1:6], replicates = rep(3,6)), nperlane=4, search.surrounding = 2, seed=1)
  expect_true(length(unique(des$Design$design$lane))==5)
  expect_true(length(unique(des$suggestedDesign$design$lane))==3)
})




