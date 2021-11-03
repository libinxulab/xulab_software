## CCSbase

This is the code for building the combined CCS database and corresponding CCS prediction model hosted at 
[CCSbase.net](http://www.ccsbase.net). Previous versions of this database and predictive model were built using code 
from https://github.com/dylanhross/c3sdb, but future database and predictive model builds will be performed using this 
repository.

<hr>

### build_c3sdb
This module contains a set of scripts that are used for building the combined CCS database (`C3S.db`) from cleaned
source datasets (housed in `build_c3sdb/cleaned_data`) in JSON format. `build_c3sdb/__main__.py` is the primary script
that coordinates the high-level steps necessary to build the database:
1. Initialize `C3S.db` with empty tables (`master`, `mqns`, and `predicted`)
2. Fill the `master` table with data from the cleaned source datasets
3. Try to fill in any missing SMILES structures by searching PubChem/LipidMAPS, or through structure generation
4. Compute molecular descriptors (MQNs) for all entries with SMILES structures, add to `mqns` table
5. Annotate entries with rough chemical classification labels (lipid, peptide, carbohydrate, or small molecule)


The source datasets to include in the database can be changed by editing a line in the `main` function of 
`build_c3sdb/fill_db_from_src.py`:
```python
# source datasets (use all defined in reference_info.py)
dsets = [_['src_tag'] for _ in ref_info]
```
_The default behavior is to include all of the source datasets that are defined in `build_c3sdb/reference_info.py`. To
change this, edit the above line to reflect the explicit list of src_tags to include in the database._


#### Usage
The database build scripts require external libraries `rdkit` and `requests` to be installed. `rdkit` can be installed
following instructions [here](https://www.rdkit.org/docs/Install.html), and `requests` is installable via Python's
built-in package manager, `pip`. The build process is performed by calling the module directly from the command line:
```bash
python3 -m build_c3sdb | tee build.log
```
As the build process proceeds, information about the build process is printed to the console so it is often a good idea
to record this information for future reference (_e.g._ using `tee` as in the example above). `C3S.db` is produced in
the current working directory, and if it already exists it will be overridden. 

<hr>

### train_prediction_model
This module contains a set of scripts that produce a CCS prediction model using the data in `C3S.db` for model 
training. The model is constructed in the same general fashion as described in our 
[original paper](https://pubs.acs.org/doi/10.1021/acs.analchem.9b05772), using K-means clustering for untargeted 
classification of compounds into structurally related groups, then training separate SVR models on each of the 
individual clustered data sets. 

#### Usage
Prediction model training reqiures the external libraries `numpy`, `matplotlib`, and `scikit-learn` to be installed.
All of these libraries are installable via Python's built-in package manager, `pip`. Model training is performed by
calling this module directly, with the path to the combined CCS database provided as an argument:

```bash
python3 -m train_prediction_model C3S.db
```

Model training can take a long time especially when using a lot of hyperparameter combinations, so it is best to do it 
on a computer with a lot of CPU threads. After model training has been completed, several files will be produced in the 
working directory:
* `kmcm_svr_final_metrics.json` - CCS prediction performance metrics for training/test datasets
* `kmcm_svr_final_metrics.png` - plot of the CCS prediction performance metrics for training/test datasets
* `kmcm_svr_final.pickle`, `OHEncoder.pickle`, `LEncoder.pickle`, `SScaler.pickle`  - individual instances of the
components of the prediction model, saved in Python's serialized binary format, these are all necessary for deployment
of CCS prediction model to the website 

There are some parameters that control aspects of model training which are defined in 
`train_prediction_model/config.py`:
```python
# pRNG seed for deterministic results in stochastic steps (e.g., splitting training/test set data)
seed = 1234

# number of individual clusters to split the dataset into
n_clusters = [5]

# SVR hyperparameters C and gamma
C = [100, 1000]
gamma = [0.001, 0.01]
```

#### `C3SD`: an interface for manipulating data from `C3S.db`
The `C3SD` object (defined in `train_prediction_model/C3SData/data.py`) provides an interface for manipulating data
from the combined CCS database. Briefly, the object can construct a combined dataset from specified sources, compute
features from MQNs stored in the database, perform encoding of MS adducts, center and scale feature data, and split
into training/test data sets. This object can be used for further development of CCS prediction outside of the standard
method described in the original paper. Detailed documentation is available directly from the module:

```python
from train_prediction_model.C3SData.data import C3SD
help(C3SD)
```

<hr>

### `generate_pubchem_csv.py`
This script takes the data from `C3S.db` and generates the two files that are needed for uploading the data to PubChem
(`ccsbase_pubchem.csv` and `ccsbase_pubchem_substance.csv`). 

#### Usage
```bash
python3 generate_pubchem_csv.py
```
