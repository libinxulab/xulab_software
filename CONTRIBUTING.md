## Before Contributing
Before contributing to this repository, you may either use the libinxulab credentials to log in to GitHub and add your
personal GitHub account as a contributor so you can work on changes directly in this repository, or you can fork the 
repository to your personal account and work on your changes there. As a general practice, do not make commits directly
to the `main` branch. The `main` branch is meant for "production" code, _i.e._ complete, documented, and tested code. 
When you are working on a new utility or making edits, begin by creating a new branch from main
with a descriptive name. As you work on your changes, you can make commits to this branch (using descriptive commit 
messages), and then when you are satisfied with your changes you can make a pull request to add the changes from this 
branch into the `main` branch. You can optionally ask that someone else review your changes before they are merged into 
`main`, or if there are no conflicts then you can merge your new branch into `main` then delete your new branch. 

## Repository Layout
The top level directories of this repository are organized into thematic groups of related utilities (_e.g._ `rnaseq/` 
and `imms/`), with lower level directories as needed.  If you are adding a new utility to the repository, add it to the 
appropriate top-level (and lower-level if applicable) directory, or create a new top-level directory if none of the 
existing ones is appropriate for the utility you are adding.

## Documentation
The primary documentation for this repository is in the form of README files (in markdown format) that describe how to
use the individual utilities. There are two levels of README files: top-level and utility-level. The top-level README 
file shows a simplified layout of the repository and serves as an index with links to the utility-level README files 
which describe the usage of single or groups of related utilities. An example of a utility-level README file is in
[rnaseq/analysis/README.md](rnaseq/analysis/README.md), which has documentation for both of the utilities (in this case,
R scripts) housed in that directory. In general, you should create a new README file whenever you are adding a
utility in a new top- or lower level directory and if you are adding a utility to a directory that already has a README
file you should edit that file to reflect your newly added utility. If you are only making changes to an existing 
utility, then simply edit the existing README file to reflect those changes. 

## Example 
Let's say we are going to add a couple of Python scripts that are useful for processing MALDI-IM-MS imaging data:
* `align_slices.py`
* `single_feature_abundance.py`

The first step before making any changes to the code would be to create a new branch (based on the most up to date 
version of the `main` branch), which we could name something like `add_maldi_imaging_scripts`.

These utilities belong in the `imms/` top-level directory, so we will add a new utility-level directory
(`imms/maldi_imaging/`) along with corresponding utility-level documentation (`imms/maldi_imaging/README.md`).
We can add the new scripts to the utility-level directory to give the final structure:
```
imms/
    maldi_imaging/
        README.md
        align_slices.py
        single_feature_abundance.py
``` 
Next, we need to write the appropriate utility-level documentation (`imms/maldi_imaging/README.md`) that describes the 
usage of the two scripts we have added in detail. An example layout (in markdown format) follows:
```markdown
## xulab_software/imms/maldi_imaging

### `align_slices.py`
Aligns two slice MS images

#### Requires:
* `slice1.npy`, `slice2.npy` (MS images for the two slices to align, numpy format)

#### Produces:
* plot of alignment result
* translation instructions for alignment (slice2 relative to slice1)

#### Usage:
1. ...
2. ...
3. ...

### `single_feature_abundance.py`
Computes the abundance of a single feature from a MS image

#### Requires:
* image.npy (MS image, numpy format)

#### Produces:
* plot of feature abundance in MS image

#### Usage:
1. ...
2. ...
3. ...
```

Which renders into:
<hr>

## xulab_software/imms/maldi_imaging

### `align_slices.py`
Aligns two slice MS images

#### Requires:
* `slice1.npy`, `slice2.npy` (MS images for the two slices to align, numpy format)

#### Produces:
* plot of alignment result
* translation instructions for alignment (slice2 relative to slice1)

#### Usage:
1. ...
2. ...
3. ...

### `single_feature_abundance.py`
Computes the abundance of a single feature from a MS image

#### Requires:
* image.npy (MS image, numpy format)

#### Produces:
* plot of feature abundance in MS image

#### Usage:
1. ...
2. ...
3. ...
<hr>


Next, we need to link to the utility-level documentation from the top-level `README.md` by adding an entry under the
`imms/` section as follows:
```markdown
* __imms/__
  * [dhrmasslynxapi](imms/dhrmasslynxapi/README.md)
  * [ccsbase](imms/ccsbase/README.md)
  * [maldi_imaging](imms/maldi_imaging/README.md)
    * `align_slices.py`
    * `single_feature_abundance.py`
```
*In this example we also list the individual scripts below the link to the utility-level documentation because this 
utility is comprised of two standalone scripts, rather than a larger or more complex module like `dhrmasslynxapi`.*

Finally, we can commit all of the changes to our `add_maldi_imaging_scripts` branch, and when we are satisfied with the 
state of the code and documentation we can create a pull request to add those changes into the `main` branch. 
 