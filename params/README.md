# Parameter Generation

## Dependencies
- [`sage-math`](https://www.sagemath.org/)
- [`lattice-estimator`](https://github.com/malb/lattice-estimator) (which we ran with commit `352ddaf`) s.t. the folder structure looks like this.
  The lattice-estimator is a submodule of this repository and needs to be instantiated via `git submodule update --init` before running the script.

Folder structure:
- `benches`
- `estimator` (contains the repository of the `lattice-estimator`)
  - `docker`
  - `docs`
  - `estimator`
  - `doctest.sh`
  - ...
- `params`
  - `param_selection.ipynb`
  - `README.md`
- `src`
- ...

This script runs in a jupyter notebook, which can be started with
```bash
sage -n jupyter
```
if `sage` is available (if it's installed in a `conda` environment, then the appropriate environment needs to be activated).
Inside the jupyter notebook, the cells need to be executed in order to ensure all functions are defined.
The execution of `is_secure(rough=False)` may take up to 30 minutes to compute the concrete hardness of LWE.
