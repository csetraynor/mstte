#------------------------------
# Below are code chunks taken from the 'rstanarm' R package, obtained
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

.onLoad <- function(libname, pkgname) { # nocov start
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
} # nocov end
