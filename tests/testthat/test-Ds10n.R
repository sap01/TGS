test_that('Ds10n is processed correctly', {
  testthat::skip_on_cran()

  ## Using namespace 'base' or 'pkgLoad' with 'system.file()'
  ## breaks the command
  input_dir <- system.file('extdata',
                           package = 'TGS',
                           mustWork = TRUE)

  output_dir <- base::paste(input_dir,
                            'Output_Ds10n',
                            sep = '/')

  TGS::LearnTgs(
    isfile = 0,
    input.data.filename = 'InSilicoSize10-Yeast1-trajectories.tsv',
    num.timepts = 21,
    true.net.filename = 'DREAM3GoldStandard_InSilicoSize10_Yeast1_TrueNet.RData',
    input.wt.data.filename = 'InSilicoSize10-Yeast1-null-mutants.tsv',
    is.discrete = FALSE,
    num.discr.levels = 2,
    discr.algo = 'discretizeData.2L.wt.l',
    mi.estimator = 'mi.pca.cmi',
    apply.aracne = FALSE,
    clr.algo = 'CLR',
    max.fanin = 14,
    allow.self.loop = FALSE,
    scoring.func = 'BIC',
    input.dirname = input_dir,
    output.dirname = output_dir,
    json.file = ''
  )

  pred_result <- base::load(
    base::paste(output_dir, 'Result.RData', sep = '/'))

  expec_result <-
    base::paste(input_dir, 'Result_Ds10n.RData', sep = '/')

  testthat::expect_known_output(pred_result, expec_result, print = TRUE)
})
