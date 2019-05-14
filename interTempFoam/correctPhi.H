CorrectPhi
(
    U,
    phi,
    p_rgh,
//ADDED-------->
//	T,
//<---------ADDED
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"
