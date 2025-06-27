# mapi (version 1.1.3)
# mapi (version 1.1.1)
## NEW FEATURES
2025-06-25: Faster implementation of adaptive grid resolution `MAPI_Varicell (cartesian grids only!)`
2024-02-20: Added experimental option ignore.weights in MAPI_RunOnGrid
2022-02-01: New code for worldwide dataset based on s2 spherical library
2021-08-31: Introducing adaptive grid resolution `MAPI_Varicell`
2019-10-21: New function `MAPI_Plot2` based on ggplot2; function `MAPI_Plot` now marked as deprecated.
## BUG FIXES
2025-06-25: Corrected matrix detection in `MAPI_CheckData`
2022-01-18: MAPI_CheckData does not anymore lowercase colnames as this avoided the use of an errRad column.
2021-11-24: The column 'permuts' generated some troubles when exporting results, it is now discarded from MAPI outputs.
2021-10-29: The number of ellipses (nb_ell) is corrected and exclude non-existent metric.
2021-08-31: Bug removed from `MAPI_RunOnGrid` in crs parsing
2020-04-09: Bug removed from `MAPI_RunOnGrid` for very large dataset while printing out the number of matches (thanks to Simon Dellicour for uncovering it)
