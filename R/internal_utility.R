.validate_var = \(data, target, source, detrend = FALSE) {
  coords = NULL
  if (detrend) {
    if (inherits(data, "sf")) {
      coords = as.data.frame(sdsfun::sf_coordinates(data))
    } else if (inherits(data, "SpatRaster")) {
      coords = as.data.frame(
        terra::rowColFromCell(data, seq_len(terra::ncell(data)))
      )
    }
    names(coords) = c("x", "y")
  }

  var_indices = c(abs(target[1]), abs(source[1]))
  if (inherits(data, "sf")) {
    data = sf::st_drop_geometry(data[, var_indices, drop = FALSE])
  } else if (inherits(data, "SpatRaster")) {
    data = terra::as.data.frame(data[[var_indices]], xy = FALSE, na.rm = FALSE)
  } else {
    data = data[, var_indices, drop = FALSE]
  }
  names(data) = c("target", "source")

  mat = as.matrix(data)
  if (!(typeof(mat) %in% c("integer", "double"))) {
    stop("Non-numeric values detected in input data. 
          Please remove columns such as dates, characters, or factors.", call. = FALSE)
  }
  
  res = vector("list", 2)
  res[[1]] = mat[, 1, drop = TRUE]
  res[[2]] = mat[, 2, drop = TRUE]
  
  return(res)
}
