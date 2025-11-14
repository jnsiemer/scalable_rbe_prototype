# Parameter Generation

## Dependencies
- [`sage-math`](https://www.sagemath.org/)
- [`lattice-estimator`](https://github.com/malb/lattice-estimator) with commit `352ddaf` s.t. the folder structure looks like this.

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
if `sage`'s environment is activated.
How you activate your `sage` environment depends on your installation of `sage`.
