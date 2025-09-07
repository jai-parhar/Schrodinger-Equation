# Solver for the 1D and 2D Schrodinger Equation

Computational solution to the Schrodinger equation in 1D and 2D systems. Written in MATLAB. Primarily written over the course of UBC course PHYS 410. 

Provides a computational simulation to the Schrodinger equation using an iterated Crank–Nicolson (CN) approach for the 1D solver, and an alternating-direction implicit (ADI) approach for the 2D solver.

## File contents
This will explain the contents of each file at a high level. Further documentation exists within each source file, at the beginning of each function.

#### sch\_1d_cn.m
Provides a dynamic simulation for the 1D Schrodinger equation using the CN method. Contains a few common initial wavefunctions and common potentials, with the ability to add more with very light modification.

#### sch\_2d_adi.m
Provides a dynamic simulation for the 2D Schrodinger equation using the ADI method. Contains a few common initial wavefunctions and common potentials, including the double slit, with the ability to add more with very light modification.

#### ctest_1d.m
Convergence tests on the 1D CN solver to test the effectiveness of the solver, and to ensure convergence of the result with higher and higher levels of spatial and temporal discretization.

#### ctest_2d.m
Convergence tests on the 2D ADI solver to test the effectiveness of the solver, and to ensure convergence of the result with higher and higher levels of spatial and temporal discretization.

#### barrier_movie.m
Code to generate an animated simulation of a wavepacket moving towards a potential well. Result can be seen in `barrier.mp4`.	

#### barrier_survey.m
Survey on the tunneling probability for a wavepacket at various barrier potentials.

#### well_movie.m
Code to generate an animated simulation of a wavepacket moving towards a potential well. Result can be seen in `well.mp4`.	

#### well_survey.m
Survey on the probability for a wavepacket to be found within a potential well for various potential well depths.

#### double\_slit_movie.m
Code to generate an animated simulation of a wavepacket moving towards a double slit. Result can be seen in `double_slit.mp4`.

#### well.mp4
Animated simulation of the 2D Schrodinger equation in a potential well generated using `well_movie.m`. Contains a 2D wavepacket moving towards a potential well. It is very interesting to note the wave patterns forming within the well, as well as the scattering off the well (at least to me)!

#### barrier.mp4
Animated simulation of the 2D Schrodinger equation in a potential barrier generated using `barrier_movie.m`. Contains a 2D wavepacket moving towards a potential barrier. It is very interesting to note a component of the wavefunction tunneling through the barrier, as well as it seeming to partially diffract around the barrier.

#### double_slit.mp4
Animated simulation of the 2D Schrodinger equation in a double slit potential generated using `double_slit_movie.m`. Contains a 2D wavepacket moving towards a double slit, modeled using a set of potential barriers with a very high potential, to supress tunneling. The diffraction pattern can be clearly seen beyond the barrier, as well as a high probability to instead reflect off the barrier.
