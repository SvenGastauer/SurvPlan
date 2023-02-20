#' geojson file to csv for strata description
#' @import geojsonR
#' @param fn path to geojson files
#' @param vrep name of the column containing the Stratum names, default: "StratumNam"
#' @param val2rep values to be replaced, all strata need to be in numeric format
#' @param valrep values that replace the values from val2rep
#' @param outfn path where the output csv should be saved
#' @export
#'
gojson2csv = function(fn, vrep = "StratumNam",valrep=c("1511","1512"), outfn){
  file_js = FROM_GeoJson(url_file_string = fn)

  ppp=do.call("rbind",lapply(1:length(file_js$features), FUN=function(i)data.frame(file_js$features[i])))

  if(length(val2rep)>0){
    for(i in 1:length(val2rep)){
      ppp[,vrep][ppp[,vrep]==val2rep[i]]=valrep[i]
    }
  }
  write.csv(x = data.frame(Lon=as.numeric(ppp$geometry.coordinates.1),
                           Lat=as.numeric(ppp$geometry.coordinates.2),
                           Stratum=as.numeric(ppp[,vrep])),
                           file=outfn)
}

