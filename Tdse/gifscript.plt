#!/usr/local/bin/gnuplot -persist
# set terminal postscript landscape noenhanced monochrome \
#              dashed defaultplex "Helvetica" 14
# set output 'output.ps'
# run load "gifscript.plt"
set terminal gif animate delay 10
set output 'fdanimate.gif'
do for [i=0:49]{
file = "psifd".i; 
set xlabel "----x-grid points ------>"
set ylabel "-----|psi|^2------->"
set yrange [0:7]
set title "Time Dependent Schrodinger Equation [Finite Difference Method]"
plot file u 1:2 w l,file u 1:3 w l;
}

set terminal gif animate delay 10
set output 'ftanimate.gif'
do for [i=0:49]{
file = "psift".i; 
set xlabel "----x-grid points ------>"
set ylabel "-----|psi|^2------->"
set yrange [0:7]
set title "Time Dependent Schrodinger Equation [Fourier Transformation Method]"
plot file u 1:2 w l,file u 1:3 w l;
}
#    EOF
