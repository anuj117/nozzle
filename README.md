# nozzle
this project contain two main,programx.
you have to call MAIN.m 
  inputs-(Me,n,gamma,D)
                        Me-exit mach number,n-number of characteristics , gama-specific heat capacity,D(in mm)-throat radius 
  outputs-(x,y,AR,max_t,yp)
                        x-(x-coordinate),y-(y-coordinates),AR-area ratio,max_t-max wall angle,yp-radius of exit
programx calculates prandtl meyer angle,mach number,reimann invariant for each right and left running characteristics using compatibility equations described in the paper
MAIN.m function uses output from programx.m ,using them in straight line-slope formulae to plot all expansion waves.
finally writing all obatined wall points to points.dat file .which will be generated after successful run of the program
