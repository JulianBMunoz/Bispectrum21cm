# Bispectrum21cm
C++ code for dealing with secondaries in 21-cm bispectra as in 1506.04152 (J.B. Munoz, Y.Ali-Haimoud, and M. Kamionkowski)


The file intdiscreteGLo.cpp calculates the secondary bispectra (21 shapes),
and the local shape for primordial non-gaussianity. Change BispNG to desired shape.
--z-indepdendent.

The file readdiscreteGLo.cpp reads the integrated bispectra and performs a cosmic-variance-limited
Fisher-matrix analysis, and computes the Covariance matrix (F^(-1)), to compute the minimum fNL measurable (as the amplitude of BispNG).
--sums over z.

The rest of files are support.
