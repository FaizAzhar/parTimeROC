test_that("Checking utils", {
  expect_equal(set_baseline("exponential"),1)
  expect_error(set_baseline(1), "Incorrect baseline")
  expect_equal(set_baseline("abcde"), NULL)
})
