# Wet Canopy Evaporation
Wet canopy evaporation is the removal of water from tree canopies during precipitation. The simulations here aim to model this phenomena to assess moisture fluxes in 2 dimensions.

- Cases have been developed with OpenFoam v2212 (Code is compatible with v2006 but all references to 'LAD' and 'Cd' will need changing to 'leafAreaDensity' and 'plantCd' respectively).
- All cases here use RANS turbulence modelling with a k-epsilon closure
- An incompressible single phase solver has been used: buoyantBoussinesqSimpleFoam and buoyantBoussinesqPimpleFoam
- Water is treated as a non reactive scalar tracer
- Simulations assume an infinite supply of water in a fully saturated canopy
- Moisture is generated as a source activated by the presence of a leaf area density
- Buoyancy effects on the moisture are not modelled, turbulent transport is considered the dominant effect

Convergence issues occur if a coarse leaf area density field is defined. A python code has been included in /Utilities/ for creating a file in the required format to be placed in /setups.orig/Canopy/system/setFieldsDict.

Two files, "Residuals" and "EPMGen" will be created at runtime in the results directory. These can be executed with gnuplot for monitoring residuals and the total moisture generation.




# Steady State
## Inlet/Outlet (2D)
The inlet/outlet cases use an analytical profile for the flow in the atmospheric boundary layer. The approach in this code looks at the flux from an isolated woodland patch.

## Cyclic (1D)
The cyclic cases attempt to find a steady state solution assuming an infinite canopy in the streamwise direction. The code can be adapted to create a infinite recurrence of patches and clearings to establish the impact of woodland edges.

# Transient (Not finalised)
The transient simulation uses the profiles from a steady state simulation in the absence of a woodland to initialise the flow for a transient solver.

