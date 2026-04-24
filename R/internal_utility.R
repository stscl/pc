.validate_var = \(data, target, source) {
  if (inherits(data, "sf")) {
    mat = as.matrix(sf::st_drop_geometry(data))
  } else if (inherits(data, "SpatRaster")) {
    mat = terra::values(data, mat = TRUE)
  } else {
    mat = as.matrix(data)
  }

  res = vector("list", 4)
  res[1] = mat[, target, drop = TRUE]
  res[2] = mat[, source, drop = TRUE]
  
  return(res)
}
