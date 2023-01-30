<About>

This is the program to calculate

・superconducting order parameter, -U<cc>, at every sites

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



<Referrence>

[1] Y. Nagai. J. Phys. Soc. Jpn. 89, 074703 (2020).
[2] R. Ghadimi, T. Sugimoto, K. Tanaka, and T. Tohyama. (2021). Phys. Rev. B. 104, 144511.
[3] https://qiita.com/cometscome_phys/items/ca1b36297a29433c7ff6 (about the module)