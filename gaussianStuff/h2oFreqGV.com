%chk=h2oFreqGV.chk
# opt=(tight,z-matrix) freq mp2/aug-cc-pvdz 

WaterOptGV

0 1
 O              
 H                  1            B1
 H                  1            B1    2            A1

   B1             0.96000000
   A1           104.50000006


