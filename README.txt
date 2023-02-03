<About>

This is the program to calculate

・superconducting order parameter, -U<c_↓ c_↑>, at every sites

and

・particle number, <c^dagger c>, at every sites and spins

self-consistently, using "LK-BdG" method.[1]

The Hamiltonian is of 2-dimentional superconductor when adding magnetic field h_z, and Rashba spin-orbit interaction.[2]


You can execute the main program, LKBdG.f90, by using the included CSRmod.f90 module.[3]



<Hop file>

You have to prepare the hopping indices and the unit vector in the direction.

The included hop.txt is the example. When the electron hopps from site j to site i, and the direction are x_ij, and y_ij, you write the file as follows:

i    j    x_ij    y_ij


The example file is for one of Penrose approximants.


<Added>
I found that this self-consistent calculation has an initial value dependence with respect to the pair potentials.
Therefore, it is necessary to calculate the free energy of the system
after all the eigenenergies are also calculated and compare them to determine the appropriate initial value.
Free_energy.f90 is the program to calculate it at absolute zero.

Empirically, it seems safe to assume that the value of the imaginary part corresponds to the filling.
The imaginary part is 0 at half-filling, and as the filling decreases from there, the imaginary part becomes positive,
and as the filling increases, the imaginary part becomes negative (although I have not been able to prove this mathematically).

(February 3, 2023)


<Referrence>

[1] Y. Nagai. J. Phys. Soc. Jpn. 89, 074703 (2020).
[2] R. Ghadimi, T. Sugimoto, K. Tanaka, and T. Tohyama. (2021). Phys. Rev. B. 104, 144511.
[3] https://qiita.com/cometscome_phys/items/ca1b36297a29433c7ff6 (about the module)
