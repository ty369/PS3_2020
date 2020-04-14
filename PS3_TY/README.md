## PS3
## 1a. S matrix is shown in the excel file, SMAtrix

## 1b. The PS3_balancecheck.m matlab file gives the balanced elements for the v1-v5rervse. We don't account for the bs because we dont know the amount across the bounds.

## 1c.
Vmax=kcat*[E]
Uses Km [mM]=0.0051 for testoerone, Mus musculus;
Find some substrate concentration bounds from Park et al's paper.

### The optimal flux that I calculated was approximately: 1.27 mmol/gDW-hr.


### Requirements
The ``Solve.jl`` solution script requires the ``GLPK`` package to the FBA problem. See [GLPK](https://github.com/JuliaOpt/GLPK.jl) for details.
