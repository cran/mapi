# mapi 1.0.3

* Added a `NEWS.md` file to track changes to the package.
# mapi (version 1.0.3)
## BUG FIX: The number of ellipses (nb_ell) is corrected and exclude non-existent metric.
# mapi (version 1.0.2)
## BUG FIXES
1. New function `MAPI_Plot2` based on ggplot2; function `MAPI_Plot` now marked as deprecated.
2. Bug removed from `MAPI_RunOnGrid` for very large dataset while printing out the number of matches (thanks to Simon Dellicour for uncovering it)
3. Bug removed from `MAPI_RunOnGrid` in crs parsing
