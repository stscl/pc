.validate_var = \(data, target, source) {
  if (inherits(data, "sf")) {
    mat = as.matrix(sf::st_drop_geometry(data))
  } else if (inherits(data, "SpatRaster")) {
    mat = terra::values(data, mat = TRUE)
  } else {
    mat = as.matrix(data)
  }

  if (!(typeof(mat) %in% c("integer", "double"))) {
    stop("Non-numeric values detected in input data. 
          Please remove columns such as dates, characters, or factors.", call. = FALSE)
  }

  res = vector("list", 2)
  res[[1]] = mat[, target, drop = TRUE]
  res[[2]] = mat[, source, drop = TRUE]
  
  return(res)
}
