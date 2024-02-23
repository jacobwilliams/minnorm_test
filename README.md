# minnorm_test
Example large m x n matrix (m&lt;n) for testing minimum norm linear solvers

### notes

create a conda env:

```
conda create -c conda-forge --prefix ./env fpm
conda activate ./env
```

set up (windows) env:

```
run.bat
```

run a test:

```
fpm build --compiler ifort --flag "-fpp"
```