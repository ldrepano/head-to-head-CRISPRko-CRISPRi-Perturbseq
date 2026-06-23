version 1.0 
workflow run_SCEPTRE{
    input {
        File response_matrix
        File grna_matrix
        File grna_target_data_frame
    }
    call run_SCEPTRE_singleton_lowmoi{
        input:
            response_matrix=response_matrix,
            grna_matrix=grna_matrix,
            grna_target_data_frame=grna_target_data_frame

    }
}

task run_SCEPTRE_singleton_lowmoi{
    input {
        File response_matrix
        File grna_matrix
        File grna_target_data_frame

    }
    command<<<


    R <<RSCRIPT

    devtools::install_github("futureverse/parallelly")
    library(parallelly)
    
    print(packageVersion("parallelly"))
    print("fraction" %in% names(formals(parallelly::availableCores)))

    install.packages("devtools")
    devtools::install_github("katsevich-lab/sceptre")
    library(sceptre)
    

    print(packageVersion("parallelly"))
    print("fraction" %in% names(formals(parallelly::availableCores)))
    print(sessionInfo())

    response_matrix=as.matrix(read.table("~{response_matrix}"))
    grna_matrix=as.matrix(read.table("~{grna_matrix}"))
    grna_target_data_frame=read.table("~{grna_target_data_frame}",header=TRUE)

    print(ncol(response_matrix))
    print(ncol(grna_matrix))

    sceptre <- import_data(
      response_matrix = response_matrix,
      grna_matrix= grna_matrix,
      grna_target_data_frame = grna_target_data_frame,
      response_names= rownames(response_matrix),
      moi = "low")
    #singleton gRNA strategy will yield DE results for each guide
    sceptre<- set_analysis_parameters(
      sceptre_object = sceptre,
      discovery_pairs = construct_trans_pairs(sceptre),
      positive_control_pairs = construct_positive_control_pairs(sceptre),
      grna_integration_strategy = "singleton")
    #using the default "maximum" guide assignment method defaults to the supplied CLEANSER assignments
    sceptre<-assign_grnas(sceptre, min_grna_n_umis_threshold=1)
    #run QC, do not remove transcripts based on low counts in perturbed cells (may reflect successful perturbation)
    sceptre <- run_qc(sceptre_object = sceptre,n_nonzero_trt_thresh=0)
    # calibration check for FDR
    sceptre <- run_calibration_check(sceptre_object = sceptre,parallel = TRUE)
    #run power check: pos cons vs neg cons 
    sceptre <- run_power_check(
      sceptre_object = sceptre,
      parallel = TRUE)
    #run discovery analysis 
    sceptre <- run_discovery_analysis(
      sceptre_object = sceptre,
      parallel = TRUE)

    discovery_result <- get_result(
      sceptre_object = sceptre,
      analysis = "run_discovery_analysis")

    write.csv(discovery_result,"SCEPTRE_results.csv")

    RSCRIPT

    >>>

    runtime {
        docker: "us.gcr.io/landerlab-atacseq-200218/mehta_sceptre:latest"
        maxRetries:5
        memory: "200G"
    }

    output{
        File SCEPTRE_results="SCEPTRE_results.csv"
    }

}