%nproc=20
%mem=10000MB

#p mp2/aug-cc-pvdz opt=z-matrix freq

h2o

0  1
H 
O 1 R1
H 1 R1p     2 A1
    Variables:
R1=1.0
R1p=2.0
A1=120.0