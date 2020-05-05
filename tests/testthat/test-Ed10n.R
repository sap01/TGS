test_that('Ed10n is processed correctly', {
  testthat::skip_on_cran()

  ## Using namespace 'base' or 'pkgLoad' with 'system.file()'
  ## breaks the command
  input_dir <- system.file('extdata',
                           package = 'TGS',
                           mustWork = TRUE)

  output_dir <- base::paste(input_dir,
                            'Output_Ed10n',
                            sep = '/')

  TGS::LearnTgs(
    isfile = 0,
    json.file = '',
    input.dirname = input_dir,
    input.data.filename = 'edi-data-10n.tsv',
    num.timepts = 21,
    true.net.filename = 'edi.net.10.adj.mx.RData',
    input.wt.data.filename = '',
    is.discrete = FALSE,
    num.discr.levels = 2,
    discr.algo = 'discretizeData.2L.Tesla',
    mi.estimator = 'mi.pca.cmi',
    apply.aracne = FALSE,
    clr.algo = 'CLR',
    max.fanin = 14,
    allow.self.loop = TRUE,
    scoring.func = 'BIC',
    output.dirname = output_dir
  )

  pred_result <- base::load(
    base::paste(output_dir, 'Result.RData', sep = '/'))

  expec_result <-
    base::paste(input_dir, 'Result_Ed10n.RData', sep = '/')

  testthat::expect_known_output(pred_result, expec_result, print = TRUE)
})
