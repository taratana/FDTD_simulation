reset
set terminal pngcairo font "Times-New-Roman,18" enhanced 
set output "test1.png"
set datafile separator ','

set pm3d map

set size ratio 1
set xrange [0:1.5]
set yrange [0:1.5]

splot "E-n=050,t=1.24e-03.csv" using ($1*1e-6):($2*1e-6):3 with pm3d notitle

###########################
reset
set terminal pngcairo font "Times-New-Roman,18" enhanced 
set output "test2.png"
set datafile separator ','

set pm3d map

set size ratio 1
set xrange [0:1.5]
set yrange [0:1.5]

splot "E-n=100,t=2.43e-03.csv" using ($1*1e-6):($2*1e-6):3 with pm3d notitle
###########################
reset
set terminal pngcairo font "Times-New-Roman,18" enhanced 
set output "test3.png"
set datafile separator ','

set pm3d map

set size ratio 1
set xrange [0:1.5]
set yrange [0:1.5]

splot "E-n=150,t=3.62e-03.csv" using ($1*1e-6):($2*1e-6):3 with pm3d notitle
###########################
reset
set terminal pngcairo font "Times-New-Roman,18" enhanced 
set output "test4.png"
set datafile separator ','

set pm3d map

set size ratio 1
set xrange [0:1.5]
set yrange [0:1.5]

splot "E-n=200,t=4.81e-03.csv" using ($1*1e-6):($2*1e-6):3 with pm3d notitle
###########################
reset
set terminal pngcairo font "Times-New-Roman,18" enhanced 
set output "test5.png"
set datafile separator ','

set pm3d map

set size ratio 1
set xrange [0:1.5]
set yrange [0:1.5]

splot "E-n=250,t=6.01e-03.csv" using ($1*1e-6):($2*1e-6):3 with pm3d notitle