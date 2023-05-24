# Wet Canopy Evaporation
Wet canopy evaporation is the removal of water from tree canopies during precipitation. The simulations here aim to model this phenomena to assess moisture fluxes in 2 dimensions.

- Cases have been developed with OpenFoam v2212
- All cases here use RANS turbulence modelling with a k-epsilon closure
- An incompressible single phase solver has been used: buoyantBoussinesqSimpleFoam and buoyantBoussinesqPimpleFoam
- Water is treated as a non reactive scalar tracer
- Simulations assume an infinite supply of water in a fully saturated canopy
- Moisture is generated as a source activated by the presence of a leaf area density
- Buoyancy effects on the moisture are not modelled, turbulent transport is considered the dominant effect


# Steady State
## Cyclic (Not finalised)
The cyclic cases attempt to find a steady state solution assuming an infinite canopy in the streamwise direction. The code can be adapted to create a infinite recurrence of patches and clearings to establish the impact of woodland edges.

## Inlet/Outlet (Not finalised)
The inlet/outlet cases use an analytical profile for the flow in the atmospheric boundary layer. The approach in this code looks at the flux from an isolated woodland patch.

# Transient (Not finalised)
The transient simulation uses the profiles from a steady state simulation in the absence of a woodland to initialise the flow for a transient solver.

