process SCRANK {
  """
  Runs an Rscript on Seurat Object
  """

  label 'process_medium'
  label "r_scrank"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://diegomscoelho/scrank:latest':
            'docker.io/diegomscoelho/scrank:latest' }"

  input:
    path obj
    val column
    val species
    val target

  output:
    path "*_${species}_${target}*.rds"

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    Rscript ${workflow.projectDir}/bin/scrank.R ${obj} ${column} ${species} ${target} 
    """
}
