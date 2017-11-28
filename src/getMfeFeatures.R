# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

getMfeFeatures = function(data) {

  common.summary = c("kurtosis", "max", "mean", "median", "min", "sd", "skewness", "var", "hist")

  cat("   - mfe general features \n")
  general = tryCatch({
    unlist(mfe::mf.general(formula = as.formula("Class ~ ."), data = data, 
      features = "all", summary = common.summary))
  }, error = function(err) {
    cat("    * got some error - returning empty vector ... \n")
    print(err)
    return(numeric(0))
  })

  cat("   - mfe statistical features \n")
  statistical = tryCatch({
    unlist(mfe::mf.statistical(formula = as.formula("Class ~ ."), data = data, 
      features = "all", summary = common.summary))
 }, error = function(err) {
    cat("    * got some error - returning empty vector ... \n")
    print(err)
    return(numeric(0))
  })

  cat("   - mfe model based features \n")
  model.based = tryCatch({
    unlist(mfe::mf.model.based(formula = as.formula("Class ~ ."), data = data, 
      features = "all", summary = common.summary))
  }, error = function(err) {
    cat("    * got some error - returning empty vector ... \n")
    print(err)
    return(numeric(0))
  })

  cat("   - mfe info theo features \n")
  infotheo = tryCatch({
    unlist(mfe::mf.infotheo(formula = as.formula("Class ~ ."), data = data, 
      features = "all", summary = common.summary))
  }, error = function(err) {
    cat("    * got some error - returning empty vector ... \n")
    print(err)
    return(numeric(0))
  })

  cat("   - mfe discriminant features \n")
  discriminant = tryCatch({
    unlist(mfe::mf.discriminant(formula = as.formula("Class ~ ."), data = data, 
      features = mfe::ls.discriminant()[-8], summary = common.summary)) 
      # sdration raises a segmentation fault on server
 }, error = function(err) {
    cat("    * got some error - returning empty vector ... \n")
    print(err)
    return(numeric(0))
  })

  cat("   - mfe landmarking features \n")
  landmarking = tryCatch({
    unlist(mfe::mf.landmarking(formula = as.formula("Class ~ ."), data = data, 
      features = "all", summary = common.summary))
  }, error = function(err) {
    cat("    * got some error - returning empty vector ... \n")
    print(err)
    return(numeric(0))
  })

  # final output
  obj = list(general = general, statistical = statistical, model.based = model.based, 
    infotheo = infotheo, discriminant = discriminant, landmarking = landmarking)

  return(obj)
}

dataName = list.files(path = "../datasets/pca/", pattern = '*.csv', recursive=T)
dataPath = paste('../datasets/pca/', dataName, sep='')

result = list()
i=0
for(f in dataPath){
    i= i+1
    print(i)
    data = read.csv2(f)
    obj = getMfeFeatures(data)

    save(obj,file=paste('../output/meta-feature/',
        paste(strsplit(dataName[i],'.csv')[[1]],
              '.meta-features.RData',sep='')
              ,sep=''))

}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
