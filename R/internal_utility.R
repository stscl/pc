.validate_var = \(data, target, source = NULL, detrend = FALSE) {
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
  
  var_indices = abs(target[1])
  if (!is.null(source)) var_indices = c(var_indices, abs(source[1]))
  
  if (inherits(data, "sf")) {
    data = sf::st_drop_geometry(data[, var_indices, drop = FALSE])
  } else if (inherits(data, "SpatRaster")) {
    data = terra::as.data.frame(data[[var_indices]], xy = FALSE, na.rm = FALSE)
  } else {
    data = data[, var_indices, drop = FALSE]
  }

  if (is.null(source)) {
    names(data) = "target"
  } else {
    names(data) = c("target", "source")
  }

  if (!all(apply(data, 2, typeof) %in% c("integer", "double"))) {
    stop("Non-numeric values detected in input data. 
          Please remove columns such as dates, characters, or factors.", call. = FALSE)
  }

  res = vector("list", 2)
  if (is.null(coords)) {
    res[[1]] = data[, 1, drop = TRUE]
    if (!is.null(source)) res[[2]] = data[, 2, drop = TRUE]
  } else {
    data = as.data.frame(cbind(data, coords))
    res[[1]] = sdsfun::rm_lineartrend("target~x+y", data = data)
    if (!is.null(source)) res[[2]] = sdsfun::rm_lineartrend("source~x+y", data = data)
  }

  return(res)
}
