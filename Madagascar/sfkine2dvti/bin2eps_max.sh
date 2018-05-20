#!/bin/sh

n1=400
n2=400

d1=5
d2=5


f='e0.5d0.5'

 inp='p'$f'.dat'
outp='p'$f'.eps'


label2='width[m]'
label1='depth[m]'

title=' '

#red blue green
cchar='bclip=2 wclip=0 brgb=1,1.0,1.0 grgb=0,0.0,0.0 wrgb=1.0,1.0,1.0'

#white red blue
#cchar='wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0'

#migration
psimage < $inp n1=$n1 d1=$d1 n2=$n2 d2=$d2 \
       label1=$label1 label2=$label2 labelsize=24  \
       width=8 height=8 $cchar  perc=99  \
      > $outp



gimp $outp &
