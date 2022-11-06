# MORBID
 Calculation of rovibrational energies for a triatomic molcule in an isolated electronic state


Explanations:

The calculations are done for one J value at a time.
The effects of the non-zero electron spin in non-singlet
electronic states are neglected so in reality, J=N.
Search for "ZERO POINT" to find

0************  ZERO POINT ENERGY  ************

            E0 =   2711.21736


0************  J =   0 ENERGIES  ************


to find the start of the energy tables (the zero
point energy is given relative to the potential
minimum) and for "J =" to find the start of the tables
for subsequent J values.

For each J value, there is a column for each symmetry
(A',A'') in the C_s(M) molecular symmetry group.
The energies with a given J value and a given symmetry
are labeled with the "non-good" quantum numbers:

KA  V2  NS   

These quantum numbers are obtained from the basis function
with the largest contribution to the eigenfunction in
question. KA is the usual $K_a$, V2 is $v_2$, the principal
bending quantum number appropriate for a triatomic molecule
with a bent equilibrium structure (see Section 17.5.2 of
P.R. Bunker and P. Jensen, Molecular Symmetry and 
Spectroscopy, NRC Research Press, Ottawa, 1998). The NS
refers to the following table in the output file:

0     ************ A1 STRETCHING FUNCTIONS ************

         FCT. #        N1        N3            ENERGY/CM-1


              1         0         0              2565.89331
              2         0         1              3836.49142
              3         0         2              5090.10808
              4         1         0              6281.45001
              5         0         3              6326.91951
              6         0         4              7544.57681
              7         1         1              7554.14147
              8         0         5              8749.22967
              9         1         2              8806.09251
             10         2         0              9821.21858
             11         0         6              9935.75386
             12         1         3             10042.64425


FCT # = NS+1 (sorry about the +1).
For a given NS value, the table now gives the N1, N3 values
for the Morse oscillator product state with the largest 
contribution to the eigenfunction. The ordinary harmonic 
oscillator quantum numbers v1 and v3 can be obtained from
N1+N3 = v1+v3.

The K_c quantum number can be obtained from symmetry 
analysis (see Chapter 12 of P.R. Bunker and P. Jensen, 
Molecular Symmetry and Spectroscopy, NRC Research Press, 
Ottawa, 1998) in conjunction with the fact that K_a + K_c = 
J or J+1. I believe that K_c is even for A' symmetry and 
odd for A'' symmetry but that should be checked.


