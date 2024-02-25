Example large `m` x `n` matrix (`m < n`) for testing minimum norm linear solvers

Require FPM.

Example: `fpm run --profile release --flag "-lblas -llapack"`



### Windows notes

create a conda env:

```
conda create -c conda-forge --prefix ./env fpm
conda activate ./env
```

set up (windows) env:

```
setup.bat
```

run a test:

```
fpm run --compiler ifort --profile release --flag "/fpp /Qmkl:sequential"
```