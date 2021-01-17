from MDAnalysis import Universe,Merge,Writer
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSD,RMSF
import pandas as pd
from matplotlib.pyplot import plot,show,xlabel,ylabel,legend,suptitle

def preliminaryAnalysis(obj, top):
    aligner1 = align.AlignTraj(obj,obj,verbose=True,in_memory=True)
    cAlpha = obj.select_atoms("name CA")

    #computing Radius of gyration
    print("Computing Radius of gyration along the trajectory:")
    gyro_list = list()
    for ts in obj.trajectory:
        gyro_list.append(obj.atoms.radius_of_gyration())
    s_gyro = pd.Series(gyro_list)

    #computing rmsd with first frame as reference
    print("Computing c-alphas RMSD with the first frame as reference:")
    rmsd1 = RMSD(cAlpha,verbose=True).run()
    rmsd_df = pd.DataFrame(rmsd1.rmsd)

    #computind rmsf
    print("Computing c-alphas RMSF:")
    ref_coordinates = obj.trajectory.timeseries(asel=cAlpha).mean(axis=1)
    ref = Merge(cAlpha).load_new(ref_coordinates[:,None, :], order="afc")
    re_align = align.AlignTraj(obj,ref,select="name CA").run()
    # need to write the trj to disk (issue 15)?
    with Writer("rmsfit.xtc", n_atoms=obj.atoms.n_atoms) as w:
        for ts in obj.trajectory:
            w.write(obj.atoms)
    #creating the fitted trj
    rmsfObj = Universe("rmsfit.xtc",top)
    #backboneRef = rmsfObj.select_atoms("backbone")
    rmsf = RMSF(cAlpha,verbose = True).run()
    rmsf_df = pd.DataFrame(rmsf.rmsf,index=cAlpha.resnums)

    return s_gyro,rmsd_df,rmsf_df


if __name__ == "__main__":
    obj1 = Universe("MUTATED.pdb","final_whole_noPBC.xtc")
    obj2 = Universe("open_V_ps.pdb","open_450ns_protein_only.dcd")
    results1 = preliminaryAnalysis(obj1,"MUTATED.pdb")
    results2 = preliminaryAnalysis(obj2,"open_V_ps.pdb")

    #initializing groups
    compGroups = ["Wild Type","Mutated"]
    comp_gyro = pd.DataFrame(columns=compGroups)
    comp_rmsd = pd.DataFrame(columns=compGroups)
    comp_rmsf = pd.DataFrame(columns=compGroups)


    #filling dataframes and checking
    comp_gyro["Wild Type"] = results2[0]
    comp_gyro["Mutated"] = results1[0]
    comp_rmsd["Wild Type"] = results2[1][2]
    comp_rmsd["Mutated"] = results1[1][2]
    comp_rmsf["Wild Type"] = results2[2][0]
    comp_rmsf["Mutated"] = results1[2][0]
    print(comp_rmsd)

    #plotting resulting dataframes
    plot(comp_gyro["Wild Type"],"b")
    plot(comp_gyro["Mutated"],"r")
    xlabel("$[10^{-1}ns]$")                          
    ylabel("[Å]")
    legend(["Open wild type","Open mutated"])
    suptitle("Gyration radius")
    show()
    plot(comp_rmsd["Wild Type"],"b")
    plot(comp_rmsd["Mutated"],"r")
    xlabel("$[10^{-1}ns]$")
    ylabel("[Å]")
    legend(["Open wild type","Open mutated"])
    suptitle("RMSD (first frame)")
    show()
    plot(comp_rmsf["Wild Type"],"b")
    plot(comp_rmsf["Mutated"],"r")
    xlabel("Residue number")
    ylabel("[Å]")
    legend(["Open wild type","Open mutated"])
    suptitle("RMSF")
    show()