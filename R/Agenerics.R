register_generic = \(name, def = NULL) {
  if (!methods::isGeneric(name)){
    if (is.null(def)) {
      def = eval(bquote(function(data, ...) standardGeneric(.(name))))
    }
    methods::setGeneric(name, def)
  }
}

for (gen in c("pc", "ops")) {
  register_generic(gen)
}
