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

