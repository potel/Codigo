set terminal epslatex size 6.0,4.0
set output 'prop_pos_local.tex' 
set style line 1 lc rgb 'red' lw 2.0
set style line 2 lc rgb 'blue' dt 2 lw 3.0 
set xlabel 'r [fm]'
set ylabel 'S(r,r) [MeV$^{-1}$fm$^{-3}$]'
set title 'Spectral function at positive energy using only local part of DOM potential'
set label at graph 0.5,0.7 '$\ell = 1$, $j = 1.5$, $E = 20$'
set key right 
set xrange [0.5:12]
plot    "prop_pos_local.txt" w l ls 2 title 'differential equation method',\
        "prop_pos_local.txt"  u 1:3 w l ls 1 title 'k-space method'
