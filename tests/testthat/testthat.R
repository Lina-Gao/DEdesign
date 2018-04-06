library(testthat)
library(DEdesign)
library(dplyr)

context("design")

test_that("trt4 rep4 Latin square",{
  des <- gendesign(trts=letters[1:4], reps=rep(4,4), nperlane=NULL)
  expect_true(all.equal(
    des$Design$design %>% filter(lane==1) %>% arrange(adp) %>% .$trt,
    c("d", "a", "b", "c")))
  expect_true(all.equal(
    des$Design$design %>% filter(lane==2) %>% arrange(adp) %>% .$trt,
    c("a", "b", "c", "d")))
  expect_true(all.equal(
    des$Design$design %>% filter(adp==3) %>% arrange(lane) %>% .$trt,
    c("b", "c", "d", "a")))
  expect_true(all.equal(
    des$Design$design %>% filter(adp==4) %>% arrange(lane) %>% .$trt,
    c("c", "d", "a", "b")))
  expect_true(all.equal(
    des$Design$BlocksEfficiency$Efficiency %>% as.character %>% as.numeric,
    c(1,1)))

  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(lane==1) %>% arrange(adp) %>% .$trt,
    c("d", "a", "b", "c")))
  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(lane==2) %>% arrange(adp) %>% .$trt,
    c("a", "b", "c", "d")))
  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(adp==3) %>% arrange(lane) %>% .$trt,
    c("b", "c", "d", "a")))
  expect_true(all.equal(
    des$suggestedDesign$design %>% filter(adp==4) %>% arrange(lane) %>% .$trt,
    c("c", "d", "a", "b")))
  expect_true(all.equal(
    des$suggestedDesign$BlocksEfficiency$Efficiency%>% as.character %>% as.numeric,
    c(1,1)))
})
