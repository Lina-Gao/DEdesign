library(testthat)
library(DEdesign)
library(dplyr)

context("design")

test_that("treatment4 replicate4 Latin square",{
  des <- gendesign(treatments=letters[1:4], replicates=rep(4,4), nperlane=NULL, seed = 1, searches = 100)
  expect_true(all.equal(
    des$Design$design %>% filter(lane==1) %>% arrange(adapter) %>% .$treatment %>% as.character(),
    c("a", "b", "c", "d")))
  expect_true(all.equal(
    des$Design$design %>% filter(lane==2) %>% arrange(adapter) %>% .$treatment %>% as.character(),
    c("d", "c", "b", "a")))
  expect_true(all.equal(
    des$Design$design %>% filter(adapter==3) %>% arrange(lane) %>% .$treatment %>% as.character(),
    c("c", "b", "a", "d")))
  expect_true(all.equal(
    des$Design$design %>% filter(adapter==4) %>% arrange(lane) %>% .$treatment %>% as.character(),
    c("d", "a", "b", "c")))
  expect_true(all.equal(
    des$Design$BlocksEfficiency$`D-efficiency` %>% as.character %>% as.numeric,
    c(1,1)))

  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(lane==1) %>% arrange(adapter) %>% .$treatment %>% as.character(),
    c("a", "b", "c", "d")))
  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(lane==2) %>% arrange(adapter) %>% .$treatment %>% as.character(),
    c("d", "c", "b", "a")))
  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(adapter==3) %>% arrange(lane) %>% .$treatment %>% as.character(),
    c("c", "b", "a", "d")))
  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(adapter==4) %>% arrange(lane) %>% .$treatment %>% as.character(),
    c("d", "a", "b", "c")))
  expect_true(all.equal(
    des$suggestedDesign$BlocksEfficiency$`D-efficiency`%>% as.character %>% as.numeric,
    c(1,1)))
})
