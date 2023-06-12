# smoothic 1.1.0

* minor update to the `smoothic()` function to include a vector of the maximum iterations to
be performed at each epsilon-telescope step (computationally advantageous)

* addition of the `plot_effects()` plotting function to plot the model-based conditional
density curves for different covariate combinations

* addition of the `plot_paths()` plotting function to plot the standardized coefficient values
through the epsilon-telescope

* addition of the `predict.smoothic()`

# smoothic 1.0.0

* major update to the `smoothic()` function to include different families of distributions

* addition of the "smooth generalized normal distribution", where an additional shape
parameter is estimated relating to the kurtosis of the error distribution (shape parameter
can also be fixed at a user-supplied value)

* new option to use `nlm()` for optimization (`optimizer = "nlm"`) or to use the manually
coded Newton-Raphson method (`optimizer = "manual"`)

* addition of the Laplace distribution, which corresponds to robust regression where the
errors are heavy-tailed

* new dataset `bostonhouseprice2`, which is a corrected version of the original `bostonhouseprice` data

* new dataset `diabetes`

# smoothic 0.1.0

* initial release

* two datasets `bostonhouseprice` and `sniffer`

* automatic variable selection using the `smoothic` function

* can choose between distributional regression (multi-parameter) with `model = "mpr"` and location-only regression (single parameter) with `model = "spr"`
