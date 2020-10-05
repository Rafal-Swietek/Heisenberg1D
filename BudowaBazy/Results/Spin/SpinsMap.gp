reset
set grid xtics
set grid ytics
set grid ztics
set encoding utf8
system('mkdir -p png')
set terminal pngcairo size 1920,1080 font 'Verdana,13  '
set view map
set ticslevel 0
set title "Spin evolution on each site"
set zlabel "S_iz" font "Verdana,14 "
set xlabel "i'th site" font "Verdana,14 "
set ylabel "time" norotate offset -1,0 font "Verdana,14 "


set key box opaque 
set pm3d map
set pm3d interpolate 0,0

unset margin
set autoscale fix

do for[i=2:16:2]{
set title "Spin evolution on each element"
set output sprintf('SpinsL%d.png',i)                           
titla = sprintf('%d_d0_spinMap.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 0.0',i)
}


set palette rgbformulae 22,13,-31