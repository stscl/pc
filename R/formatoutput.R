#' @export
#' @noRd
print.rpc_res = \(x,...){
  pc = x$xmap
  if (x$bidirectional){
    pc$direction = rep(c(paste0(x$varname[1], " -> ", x$varname[2]),
                         paste0(x$varname[2], " -> ", x$varname[1])),
                       each = (nrow(pc) / 2))
  } else {
    pc$direction = paste0(x$varname[1], " -> ", x$varname[2])
  }
  names(pc)[5] = "strength"
  print(pc[,c(1,2,5,7),drop = FALSE])
}

#' @export
#' @noRd
print.pc_ops = \(x,...){
  res = x$param
  cat(paste0("The suggested E, k, tau is ", res[1], ", ", res[2], " and ", res[3]), "\n")
}
