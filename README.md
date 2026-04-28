Core code:\
**`stxs_reweight_function.ipynb`** is used for running processor (producing coffea) and testing stxs function code.\
**`stxs_functions.py`** contains all the functions necessary for quadratic fitting and likelihood scans.\
**`vbfprocessor.py`** is the processor used in **`stxs_reweight_function.ipynb`**.\
**`fit_plots.ipynb`** is used for plotting likelihood scans and reweighted STXS bins. **`stxs_plotting.ipynb`** also contains stxs plots.\
**`frankensteincfg.py`** is my modified pythia cfg. It uses CMS_SW_13_2_9.


Note: **`stxs_reweight_fit.ipynb`** only works for two operators or less at a time for now; working on generalizing to n-dimensional fits

Typical workflow (assuming Higgs VBF root samples have already been simulated, including reweight points for desired operators):
1. Creating coffea histograms:\
   In **`stxs_reweight_function.ipynb`**, call `run_samp` (for a single local .root file) or `remote_run_samp` (for multiple remote .root files).
3. Extracting quadratic dependencies:\
   In **`stxs_reweight_function.ipynb`** or importing **`stxs_functions.py`**, make sure to modify the default reweight card path in `wc_map_dict` if relevant.\
   In **`stxs_reweight_function.ipynb`** or importing **`stxs_functions.py`**, call `stxs_reweight_function` on a chosen coffea file and desired list of operators. This outputs a list of        coefficients. Examine the function to see order of coefficients.\
   In **`stxs_reweight_function.ipynb`**, call `plot_1D` to visualize quadratic dependencies.
5. Likelihood scans:\
   In **`stxs_reweight_function.ipynb`** or importing **`stxs_functions.py`**, call`stxs_fit` on a chosen coffea file and desired list of operators. This outputs a dictionary containing         various points in WC linspace. Each point in WC space also stores a bin1 xsec, bin2 xsec, bin1 chi-squared, bin2 chi-squared, and total chi-squared.
6. Visualizing likelihoods:\
   In **`fit_plots.ipynb`**, use `values_1d_plot`, `chi2_1d_plot`, and `chi2_2D_heatmap` to visualize ∆chi^2.
8. Reweighting STXS:\
   In **`fit_plots.ipynb`**, use `predicted_stxs_plot` to view the STXS binning for the sample reweighted to the best-fit WC point.
