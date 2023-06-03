# ProjetLU3ME008
Evolution thermique instationnaire: Traitement thermique d’une plaque métallique --- Projet numérique: LU3ME008

!Formule de compilation : gfortran -fdefault-real-8 -o exe Projet_mod.f90 Projet.f90
!Formule d'exécution    : ./exe
!!!GNUPLOT 
!set title 'Evolution de la temperature dans une plaque metalique'
!set xlabel 'Position (en metre)'
!set ylabel 'Temperature (en degre Celsius)'
!set autoscale
!set grid 
!plot 'graph_30.dat' w lp, 'graph_60.dat' w lp, 'graph_90.dat' w lp, 'graph_120.dat' w lp
!plot "graph_normx.dat" w lp, "graph_normta.dat" w lp
