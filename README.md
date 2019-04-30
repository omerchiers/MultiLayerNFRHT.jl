# MultiLayerNFRHT

This package computes the near-field radiative heat transfer between 2 plane parallel semi-infinite bodies.
Each body (substrate) can be composed of an arbitrary number of homogeneous layers of finite thickness.
The method uses the equation (Polder and Van Hove) for two bulk media but where the global reflection coefficients are the ones for the multilayer. These reflection coefficients are obtained using the S-matrix method, in order to get convergence in the near-field regime.

## To-do list
  * In the mid-term, move optical properties to a separate package. (done)
  * Add rt function for anisotropic materials.
  * Long term: include rt function for corrugated substrates.

## Installation
You need Julia version 0.7 or higher.

To install on julia 0.7 or higher, enter the Pkg mode by typing ]
and then to install from github

```julia
(v1.0)> add https://github.com/omerchiers/MyPhysicalConstants.jl
(v1.0)> add https://github.com/omerchiers/OpticalProperties.jl
(v1.0)> add https://github.com/omerchiers/MultiLayerNFRHT.jl
```
