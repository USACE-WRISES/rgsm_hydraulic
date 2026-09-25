library(rsconnect)

# Refresh authorization before deploying
rsconnect::connectCloudUser()

# Write over app and data files
# Only copy files that change; should mostly just be the app file
# Overwrite app.R in bundle with latest version
file.copy("code/shiny/app.R", "deploy_bundle/app.R", overwrite = TRUE)

# If need to write over all files:
# Rebuild bundle folder fresh
# unlink("deploy_bundle", recursive = TRUE)
# dir.create("deploy_bundle", recursive = TRUE)
# dir.create("deploy_bundle/output")
# dir.create("deploy_bundle/data/flow_perm",          recursive = TRUE)
# dir.create("deploy_bundle/data/yackulic2022_data",  recursive = TRUE)
# 
# file.copy("code/shiny/app.R",                        "deploy_bundle/app.R")
# file.copy("output/lcc_low_combined_mcmc_sub.csv",    "deploy_bundle/output/")
# file.copy("output/inund_lcc_list_combined.RData",    "deploy_bundle/output/")
# file.copy("output/ee_list.RData",                    "deploy_bundle/output/")
# file.copy(list.files("data/flow_perm", full.names = TRUE), "deploy_bundle/data/flow_perm/")
# file.copy("data/yackulic2022_data/abq_gage_08330000.csv",       "deploy_bundle/data/yackulic2022_data/")
# file.copy("data/yackulic2022_data/SanAcacia_gage_08354900.csv", "deploy_bundle/data/yackulic2022_data/")
# file.copy("data/ecoval2d.csv",             "deploy_bundle/data/")
# file.copy("data/q_hab_lookup_combined.csv", "deploy_bundle/data/")
# file.copy("data/san_acacia_rest_ecovals.csv", "deploy_bundle/data/")

rsconnect::deployApp(
  appDir        = "deploy_bundle",
  appPrimaryDoc = "app.R",
  appName       = "rgsm-hydraulic",
  server        = "connect.posit.cloud"
)

