from MDAnalysis import Universe
from matplotlib.pyplot import plot,show,xlabel,ylabel,suptitle,legend

obj1 = Universe("MUTATED.pdb","final_whole_noPBC.xtc")
obj2 = Universe("open_V_ps.pdb","open_450ns_protein_only.dcd")

temp1 = obj1.