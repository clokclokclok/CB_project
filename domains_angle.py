from numpy import arccos,dot,rad2deg
from numpy.linalg import norm
from MDAnalysis import Universe
from pandas import Series
from matplotlib.pyplot import plot,show,xlabel,ylabel,legend,suptitle




def getAngle(obj,sampling=1,roundOpt = False):
    angleResults = dict()
    for ts in obj.trajectory:
        #obtainig geometrical center (center of mass)*
        dom1 = obj.select_atoms("resid 0-99").center_of_mass()
        dom2 = obj.select_atoms("resid 105-170").center_of_mass()
        dom3 = obj.select_atoms("resid 173-252").center_of_mass()
        edge1 = dom3 - dom2
        edge2 = dom1 - dom2
        angle = arccos(dot(edge1,edge2)/(norm(edge1)*norm(edge2)))
        angleDegree = rad2deg(angle)
        if round:
            angleDegree = round(angleDegree,4)
        angleResults[ts.time*sampling] = angleDegree
        print(ts.time/len(obj.trajectory)*10*sampling)
    angleSerie = Series(angleResults)
    return angleSerie

if __name__ == "__main__":
    #creating universe
    obj1 = Universe("MUTATED.pdb","final_whole_noPBC.xtc")
    obj2 = Universe("open_V_ps.pdb","open_450ns_protein_only.dcd")
    print("Analyzing first trajectory...")
    result1 = getAngle(obj1)
    print("The mean angle of Mutated trajectory is:{}\n".format(result1.mean()))
    print("Analyzing secon trajectory...")
    result2 = getAngle(obj2,sampling=2000)
    print("The mean angle of Wild type trajectory is:{}\n".format(result2.mean()))

    print(result1)

    #plotting results
    plot(result1,"r")
    plot(result2,"b")
    xlabel("[ps]")
    ylabel("[deg]")
    legend(["Open mutated","Open wild type"])
    suptitle("Angle between domains")
    show()