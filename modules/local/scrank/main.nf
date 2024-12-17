process SCRANK {
  scratch true
  
  """
  Runs an Rscript on Seurat Object
  """

  label 'process_medium'
  label "r_scrank"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://diegomscoelho/scrank:v1.0':
            'docker.io/diegomscoelho/scrank:v1.0' }"

  input:
    val(obj)
    val(column)
    val(species)
    val(target)

  output:
    path("${species}_${target}.rds")

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    Rscript scrank.R '${obj}' '${column}' '${species}' '${target}' 
    """
}
