This is the 2024 version of FORTRAN code QSSP for calculating complete synthetic seismograms of a spherical earth using the normal mode theory. In comparison with its earlier versions, this version changes the input format to provide more options of selecting source types as well as their composition.

Highlights:

(1) all-in-one code for body waves, surface waves, free oscillations, tsunami for uniform ocean, infrasound/accoustic-gravity/Lamb waves for a standard atmosphere, and static solid-earth deformation as well

(2) generating Green’s function database or simulating complete seismograms for any given kinematic source model

(3) hybrid algorithm (numerical integration for low frequency / small harmonic degrees and analytical propagator algorithm for high frequency / large harmonic degrees)

(4) complex frequency technique for supressing the time-domain aliasing problem

(5) differential filter technique for suppressing numerical phases (i.e., space-domain aliasing, see Wang and Wang (2007))

Related codes

QSSPSTATIC - Co- and post-seismic viscoelastic deformation based on a spherical visco-elastic-gravitational earth model.

QSSPCOSEIS - Co-seismic static deformation based on a spherical elastic-gravitational earth model.

SPGRN - synthetic Green's function database based on a spherical elastic-gravitational earth model.

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../SourceCode

~>gfortran -o qssp2024 *.f -O3

to get the excutable code qssp2024.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

References

Wang, R., S. Heimann, Y. Zhang, H. Wang, and T. Dahm (2017). Complete synthetic seismograms based on a spherical self-gravitating Earth model with an atmos-phere-ocean-mantle-core structure. Geophysical Journal International, doi: 10.1093/gji/ggx259.

Wang, R., and H. Wang (2007), A fast converging and anti-aliasing algorithm for Green’s functions in terms of spherical or cylindrical harmonics, Geophysical Journal International, doi: 10.1111/j.1365-246X.2007.03385.x.
