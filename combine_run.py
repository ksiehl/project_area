import os
import sys

#this is a simple script to run the combine command over all points,
#over each parameter

def run_combine(par_of_int,ccw,cwww,cb):
    frozen = 'void'
    if par_of_int == "ccw":
       frozen = "cwww,cb"
    elif par_of_int == 'cwww':
       frozen = 'ccw,cb'
    elif par_of_int == 'cb':
       frozen = 'cwww,ccw'
    else:
       print "not an acceptable parameter.\n"
       sys.exit(1)
 
    point = 'cwww='+cwww + ',ccw='+ccw + ',cb='+cb

    name = '_'+par_of_int +'_val_cwww'+cwww+'_ccw'+ccw+'_cb'+cb

    print "Keeping parameters ", frozen, " all frozen, and varying \
the parameter ", par_of_int;
    print "Operating on the point ", point
    print "Name is ", name

    command = 'combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances '+frozen+' --setPhysicsModelParameters '+point+' --minimizerStrategy 2 --cminPreScan --redefineSignalPOIs '+par_of_int+' --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n '+name

    print 'running the command:\n'+command
    #os.system(command)
   

param_list = ['ccw','cwww','cb']

ccw_list = ['-4.5','-2.25','0','2.25','4.5']
cwww_list = ['-3.6','-1.8','0','1.8','3.6']
cb_list = ['-20','-10','0','10','20']

for par in param_list:
  for ccw in ccw_list:
    for cwww in cwww_list:
      for cb in cb_list:
        run_combine(par,ccw,cwww,cb)
