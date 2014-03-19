pap_sear_3d_sort
================

Implementace Shearsortu a 3DSortu pro OpenMP.

Pouziti:
-S NxN    Razeni shearsortem v 2d mrizce (pocet threadu bude NxN)
-3 NxNxN   Razeni 3dsortem v 3d mrizce (pocet threadu bude NxNxN)
-f path/to/file   Vstup neserazenych cisel oddelenych whitespacem
-o path/to/file   vystup serazenych cisel (nemusi byt pak stdout)

./dist/Debug/GNU-Linux-x86/sorty -3 3x4x0 -f test -o out
