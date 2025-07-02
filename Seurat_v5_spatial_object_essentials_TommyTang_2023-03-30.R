#####
##### SJM R script and notes
##### 2025-07-01
#####
##### Is this vizgen resource now out-of-date or unavailable? (gsutil is giving errors)
#####

# Following https://satijalab.org/seurat/articles/spatial_vignette_2.html

data.dir <- "C:/Downloads/spatial_data"

# Download the public breast vizgen cancer data

# gsutil is Google Cloud CLI, https://cloud.google.com/storage/docs/gsutil_install#install
#
# issue this command outside R (e.g. in a suitable terminal)
# gsutil -m cp -n gs://vz-ffpe-showcase/HumanOvarianCancerPatient2Slice2/cell_by_gene.csv gs://vz-ffpe-showcase/HumanOvarianCancerPatient2Slice2/cell_metadata.csv ./spatial_data/
  