from domains_angle import getAngle
from MDAnalysis import Universe
from statsmodels.graphics.tsaplots import plot_acf
from matplotlib.pyplot import show,plot,suptitle,legend
import pandas as pd
import numpy

#open
obj1 = Universe("MUTATED.pdb","final_whole_noPBC.xtc")
obj2 = Universe("open_V_ps.pdb","open_450ns_protein_only.dcd")
#closed
obj3 = Universe("N34I_mutated.gro","md_0_1_noPBC.xtc")
obj4 = Universe("closed_II_ps.pdb","closed_450ns_protein_only.dcd")

    
angleResults1 = getAngle(obj1,roundOpt=True)
angleResults2 = getAngle(obj2,sampling=2000,roundOpt=True)

angleResults1 = angleResults1[::10].reset_index(drop=True)

#angleResults3 = getAngle(obj3,roundOpt=True)
#angleResults4 = getAngle(obj4,sampling=2000,roundOpt=True)

#angleResults3 = angleResults3[::10].reset_index(drop=True)



# Creating Autocorrelation plot 
toPlot1 = pd.Series(pd.plotting.autocorrelation_plot(angleResults1))
toPlot2 = pd.Series(pd.plotting.autocorrelation_plot(angleResults2)) 

#toPlot3 = pd.Series(pd.plotting.autocorrelation_plot(angleResults3))
#toPlot4 = pd.Series(pd.plotting.autocorrelation_plot(angleResults4)) 
  

plot(toPlot1,"r")
plot(toPlot2,"b")
suptitle("Autocorrelation of angle between domains, open conformation")
legend(("Mutated","Wild type"))
show()
"""
plot(toPlot3,"r")
plot(toPlot4,"b")
suptitle("Autocorrelation of angle between domains, closed conformation")
legend(("Mutated","Wild type"))
show()"""