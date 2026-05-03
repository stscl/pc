.validate_var = \(data, target, source, detrend = FALSE) {
  var_indices = c(abs(target[1]), abs(source[1]))
  if (inherits(data, "SpatRaster")) {
    data = data[[var_indices]]
  } else {
    data = data[, var_indices, drop = FALSE]
  }

  if (detrend) {
    coords = NULL
    if (inherits(data, "sf")) {
      coords = as.data.frame(sdsfun::sf_coordinates(data))
    } else if (inherits(data, "SpatRaster")) {
      coords = terra::rowColFromCell(data, seq_len(terra::ncell(data)))
    }
  }


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
  res[[1]] = mat[, abs(target[1]), drop = TRUE]
  res[[2]] = mat[, abs(source[1]), drop = TRUE]
  
  return(res)
}
