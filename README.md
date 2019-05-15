# interTempFoam
Implementation of temperature profile into the interFoam OpenFOAM solver

## Instructions
The solver has been modified to account for the transport of a scalar field, temperature (variable T herein). This augmented solver requires the recompilation of of three major libraries in OpenFOAM (OF). These libraries are the "immiscibleIncompressibleTwoPhaseMixture", the "incompressibleTwoPhaseMixture", and "interMixingFoam". These libraries are modified and compiled before compiling the new solver, interTempFoam. These libraries are modified under new library names, often modifying a library of name "library" to something of the form "myLibrary"; this ensures that moidification and recompilation of these new libraries does not modify the library already compiled and implemented in OF. Provided below is a rough step-by-step procedure on the modification of the OF solver interFoam into interTempFoam. Currently, there is no interTempFoam publicly available on OF5x, and this solver aims to fill that void. It should be noted that, ideally, for the casting process, this would modified into interTempPhaseChangeFoam, which accounts for phase change in the liquid to a solid phase (this is also a solver that is not available in OF).

### Copy to user directory
#### User source code
Copying the 

#### User solvers

### incompressibleTwoPhaseMixture
First, since the other libraries are dependent on it, the incompressibleTwoPhaseMixture library must first be modified. OpenFOAM consists of two types of files, in general: class files (.C) and header files (.H). Within the class file, the temperature variable must be added under the "constructors" heading by adding the line

const volSCalarField& T

below the definition of the scalar field phi. As is made obvious by the definition, this defined the temperature variable T as a scalar field, applied to the entire volume of the fluid.

Because we will be solving the temperature equation, we need to define a few variables that are presently undefined in the interFoam solver: the heat capacity (Cp) and the thermal conductivity (kappa, or k). We will use the variable kappa, since k is used to refer to *some variable* in some of OF's k-epsilon turbulence models. This is done by adding the varibale to the class file by defining it as

cp1_("cp",dimensionSet(0, 2, -2, -1, 0, 0, 0), nuModel1_->viscosityProperties()), 
...
kappa_1("kappa",dimensionSet(1, 1, -3, -1, 0, 0, 0), nuModel1_->viscosityProperties()),
...

where the pattern is analogously followed for cp_2 and kappa_2. Notice that the dimensions of OF are defined as (Mass, Length, Time, Temperature, Quantity, Current, Luminous Intensity), with SI units applied to the dimensions

Below this the temperature field is again defined as the U and phi variables are. It is then seen that *nu* is defined (the kinetmatic viscosity). Following the convention implemented by *nu*, we will define *kappa*, the thermal conductivity. Recall that kappa has only been defined for species 1 and 2 above, and that a total kappa for the liquid/gas mixture must be calculated, weighted by the phase-fraction of liquid or gas in any given volume of fluid.
