test_that("multiplication works", {

  LearnTgs(
    isfile = 0,
    input.data.filename = "DmLc3E.RData",
    num.timepts = 6,
    is.discrete = TRUE,
    num.discr.levels = 2,
    mi.estimator = "mi.pca.cmi",
    apply.aracne = FALSE,
    clr.algo = "CLR",
    max.fanin = 14,
    allow.self.loop = TRUE,
    input.dirname = "location where file is stored",
    output.dirname = "location where output needs to be stored")

  LearnTgs(isfile = 0,
           input.data.filename = "",
           num.timepts = 0,
           true.net.filename = "",
           input.wt.data.filename = "",
           is.discrete = TRUE,
           num.discr.levels = 2,
           discr.algo = "",
           mi.estimator = "mi.pca.cmi",
           apply.aracne = FALSE,
           clr.algo = "CLR",
           max.fanin = 14,
           allow.self.loop = TRUE,
           scoring.func = "BIC",
           input.dirname = "",
           output.dirname = "",
           json.file = "")

  expect_equal(2 * 2, 4)
})
