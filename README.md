# MultiLayerNFRHT

This package computes the near-field radiative heat transfer between 2 plane parallel semi-infinite bodies.
Each body (substrate) can be composed of an arbitrary number of homogeneous layers of finite thickness.
The method uses the equation (Polder and Van Hove) for two bulk media but where the global reflection coefficients are the ones for the multilayer. These reflection coefficients are obtained using the S-matrix method, in order to get convergence in the near-field regime.

## To-do list
..* Change the OptProp abstract type to Material.
..* Replace singleton types for Materials by instances of Model types.
..* Add plot recipe for representing optical properties.
..* In the mid-term, move optical properties to a separate package.
..* Add rt function for anisotropic materials.
..* Long term: include rt function for corrugated substrates. 


