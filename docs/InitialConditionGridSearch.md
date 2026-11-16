# Initial Condition Grid Search
In addition to the parameter space exploration methods included in the `InitCond` and `ModelSelector` modules, there is an optional grid search module that can be used to provide initial conditions in place of or in addition to the standard approach. This module is called after the `InitCond` and before `Levmar`. The grid search is based on the initial condition grid search (ICGS) approach described in [Bennett 2010](https://scixplorer.org/abs/2010ApJ...716.1408B/abstract), where the chi2 of each model in the grid is evaluated, but optimization is only performed on the best initial conditions on the grid. This method was designed to be effective for low mass planetary microlensing events but could be used for any model class in RTModel. 

## The `ICGS` Module
This module can be called with the method `ICGS()`

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.config_ICGS(grid_dictionary = './LSGrid.json')
rtm.Reader()
rtm.InitCond()
rtm.launch_fits('PS')
rtm.ModelSelector('PS')
rtm.ICGS()
rtm.launch_fits('LS')
```

The ICGS can also be used with the `run()` method. In either case, `config_ICGS()` must be called.

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.config_ICGS(grid_dictionary = './LSGrid.json')
rtm.run()
```



In addition to calling the `Reader` and `InitCond` modules, this code configures the initial condition grid search code with `config_ICGS()` and runs it for the LS model search in addition to the standard binary lens initial conditions produced by RTModel's template search. 

The `ICGS` module creates the following files in a directory called `GridOuputs`:
- `{modeltype}ICGSbest_models.csv`: a csv containing the inital conditions selected from the grid. The file name starts with the model category the grid was run on, e.g, `LSICGSbest_models.csv`.
- `{modeltype}ICGS.csv`: a csv containing the model parameters and chi squared for every model on the grid.
- The resulting inital conditions are also added to the RTModel initial conditions file (e.g, `LSInitConds.txt`). They can be appended to the initial conditions created by the `InitCond()` module or can overwrite them (see below)

## Defining the Grid
The ICGS works best when the parameters for a simpler model can be well defined. For example, if the PS fit returns relatively accurate parameters, the LS initial condition grid search can scan over `s`, `q`, and `alpha`, perturbatively searching for good solutions while holding the single lens parameters fixed. The user can tell the module which parameters to scan and which to hold fixed by providing it with a JSON file. A very sparse example of such a file for the LS model class is displayed below.

```
{"modelcode": "LS",
"s":[0.9, 1.00,1.1],
"q":[1.0e-06,  1.0e-05, 1.0e-04, 1.0e-03, 1.0e-02],
"u0":[],
"alpha":[-3.15 ,  0.0,  3.15],
"rho":[0.001],
"tE":[],
"t0":[]
}
```
Each key corresponds to a different model parameter. A parameter list with multiple values is scanned over in the grid search, and a list with only one value value is held fixed in the grid search. An empty list indicates that RTModel should find the parameter from the earlier RTModel fits performed for a simpler model. For example, for the LS category, setting `"tE" = []` means the code will find the best fitting PS model from the PS fitting process and use it's Einstein timescale. Additionally, the code will add the chi squared of this model to the dictionary so the chi squared improvement of the searched grid models can be computed later. If you specify every parameter in the grid dictionary, you should also provide this chi2 in the dictionary. 

Currently, support is only provided for retrieving fixed values from the PS and PX fits for LS and LX searches, respectively, but support will be added in the future for triple lens models. Additionally, the code for this parameter retrieval is written in the class [`FindFixedParameters()`](https://github.com/jonhulberg/RTModel/blob/7200091e7a1a828e01e8a64e2d141aa4d800ca49/RTModel/icgs/icgs_helpers.py#L7). Interested users can add a method for any of the model categories in  `RTModel` and recompile the code. 

## Selection of Initial Conditions
After running the grid search, the `ICGS` module uses the class [`SelectInitConds()`](https://github.com/jonhulberg/RTModel/blob/5c0d2c2aa6a532cc4b9c8d0ec9ce1e4fcd9571a7/RTModel/icgs/icgs_helpers.py#L120). This class currently has methods defined for retrieving the best initial conditions from the LS and LX grid searches. The user can define their own methods for other model categories, and in the future we will add support for more categories. For LS and LX, the 20 best fitting models are selected, but repeats of the same `q` and `s` values with different values of `alpha` are ignored so that multiple solutions can be found.

## Options

For the calculation of lightcurves and their  chi squared, the `ICGS` module uses the same options provided to `LevMar`. Users may refer to [Fitting.md](https://github.com/jonhulberg/RTModel/blob/5c0d2c2aa6a532cc4b9c8d0ec9ce1e4fcd9571a7/RTModel/icgs/icgs_helpers.py#L120).

Additional options to control the grid search can be provided to the `config_ICGS` function. These options and their default values are enumerated below.
- `grid_dictionary = None`: a python grid dictionary which can be provided instead of a path to a JSON file.
- `grid_file_path = None`: a path to the grid JSON file. The default for both grid dictionary parameters is `None`, meaning that running RTModel without calling `config_ICGS()` and providing one of the two options will skip the ICGS. This is the default behavior of RTModel.
- `Tol = 0.01`, `RelTol = 0.001`: allows the user to control the `VBMicrolensing` accuracy parameters. This is currently not implemented in `Levmar` and should be used with caution.
- `overwrite = False`: Controls whether the module overwrites the initial conditions provided by the `InitCond` module (`True`) or appends to them (`False`).
