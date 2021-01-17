from MDAnalysis import Universe
from MDAnalysis.analysis import pca,align
from cosine_content import getCosineContent,getCosineContentDistriution
from matplotlib.pyplot import plot,xlabel,ylabel,show,scatter,figure,legend,suptitle
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import plotly.express as px
import pandas as pd

#creating Universe
#obj1 = Universe("N34I_mutated.gro","md_0_1_noPBC.xtc")
obj2 = Universe("MUTATED.pdb","final_whole_noPBC.xtc")
#print("The lentgh of the trajectory consist in:{} frames.".format(len(obj1.trajectory)))

#firstly i need to align the trajectory:  ___----->not true probably
print("Aligning trajectories...")
#aligner1 = align.AlignTraj(obj1, obj1, select="name CA", in_memory=True, verbose=True).run()
aligner2 = align.AlignTraj(obj2, obj2, select="name CA", in_memory=True, verbose=True).run()
#then perform the Principal components analysis:
print("Extracting principal components...")
#pc_obj1 = pca.PCA(obj1, select="name CA", align=False, mean = None, n_components= None, verbose=True).run()
pc_obj2 = pca.PCA(obj2, select="name CA", align=False, mean = None, n_components= None, verbose=True).run()



#first trj
"""
print("Total variance of MUTATED trajectory is:{}".format(sum(pc_obj1.cumulated_variance)))
print("Plotting the cumulative variance of MUTATED principal components:")
#plotting cumulated vanriance per Principal component:
plot(pc_obj1.cumulated_variance[:10])
xlabel("Principal component")
ylabel("Cumulative variance")
show()

#plotting atom pos in pc space
print("Projecting atom position in the space of firsts 8 MUTATED principal components:")
nameCA1 = obj1.select_atoms("name CA")
transformed1 = pc_obj1.transform(nameCA1,n_components=8)
df1 = pd.DataFrame(transformed1,
                   columns=["PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"])
df1["Time(ps)"] = df1.index * obj1.trajectory.dt
df1.to_csv("Mutated_closed_pca.csv")



#second trj
"""
print("Total variance of Wild Type trajectory is:{}".format(sum(pc_obj2.cumulated_variance)))
print("Plotting the cumulative variance of Wild Type principal components:")
plot(pc_obj2.cumulated_variance[:10])
xlabel("Principal component")
ylabel("Cumulative variance")
show()

#plotting atom pos in pc space
print("Projecting atom position in the space of Wild Type firsts 8 principal components:")
nameCA2 = obj2.select_atoms("name CA")
transformed2 = pc_obj2.transform(nameCA2,n_components=8)
df2 = pd.DataFrame(transformed2,
                   columns=["PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"])
df2["Time(ps)"] = df2.index * obj2.trajectory.dt
df2.to_csv("MUTATED_pca.csv")
"""
"""
"""
# axes instance
fig = figure(figsize=(6,6))
ax = Axes3D(fig)
sc = ax.scatter(df2["PC1"],df2["PC2"],df2["PC3"], 
                c=np.arange(len(df2)), s=40, marker='o', alpha=1)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")
# legend
legend(*sc.legend_elements(), bbox_to_anchor=(1.05, 1), loc=2)
show()
"""

# axes instance
fig = plt.figure(figsize=(6,6))
ax = Axes3D(fig)


sc = ax.scatter(transformed2[:,0], transformed2[:,1], transformed2[:,2], 
                c=np.arange(transformed2.shape[0]), s=40, marker='o', alpha=1)
# legend
plt.legend(*sc.legend_elements(), bbox_to_anchor=(1.05, 1), loc=2)
plt.savefig("pca_withLegend")
show()
"""
#cosine content of MUTATED
print("Computing cosine content of firsts 3 pc of the MUTATED trajectory...")
cosineCont1_PC1 = getCosineContent(df1["PC1"])
cosineCont1_PC2 = getCosineContent(df1["PC2"])
cosineCont1_PC3 = getCosineContent(df1["PC3"])
print("The cosine content of firsts 3 principal component is:\nPC1\t{}\nPC2\t{}\nPC3\t{}".format(cosineCont1_PC1,cosineCont1_PC2,cosineCont1_PC3))
"""
"""
#cosine content of WT
print("Computing cosine content of firsts 3 pc of the WT trajectory...")
cosineCont2_PC1 = getCosineContent(df2["PC1"])
cosineCont2_PC2 = getCosineContent(df2["PC2"])
cosineCont2_PC3 = getCosineContent(df2["PC3"])
print("The cosine content of firsts 3 principal component is:\nPC1\t{}\nPC2\t{}\nPC3\t{}".format(cosineCont2_PC1,cosineCont2_PC2,cosineCont2_PC3))
"""
"""
with open("cosineContent.txt","w") as cosRes:
    cosRes.write("#Results of cosine content analysis\n")
    cosRes.write("#Mutated\n")
    cosRes.write("PC1\t"+str(cosineCont1_PC1)+"\n")
    cosRes.write("PC2\t"+str(cosineCont1_PC2)+"\n")
    cosRes.write("PC3\t"+str(cosineCont1_PC3)+"\n")
"""
"""
    cosRes.write("#Wild type\n")
    cosRes.write("PC1\t"+str(cosineCont2_PC1)+"\n")
    cosRes.write("PC2\t"+str(cosineCont2_PC2)+"\n")
    cosRes.write("PC3\t"+str(cosineCont2_PC3)+"\n")
"""
"""
print("Plotting mutated trajectory on PCA's")
fig1 = px.scatter_3d(df1, x='PC1', y='PC2', z='PC3',
              color=df1.index.values, width=900,height=800)
fig1.show()
"""
"""
print("Plotting WT trajectory on PCA's")
fig2 = px.scatter_3d(df2, x='PC1', y='PC2', z='PC3',
              color=df2.index.values, width=900,height=800)
fig2.show()
"""
#computing cosine content distribution per ts for PCs
"""
PC1_mutated_dist = getCosineContentDistriution(df1["PC1"])
PC2_mutated_dist = getCosineContentDistriution(df1["PC2"])
PC3_mutated_dist = getCosineContentDistriution(df1["PC3"])
plot(PC1_mutated_dist)
plot(PC2_mutated_dist)
plot(PC3_mutated_dist)
xlabel("ts")
ylabel("Cosine content")
legend("PC1","PC2","PC3")
suptitle("Open mutated cosine content of firsts 3 principal components")
show()
"""
"""
PC1_wt_dist = getCosineContentDistriution(df2["PC1"])
PC2_wt_dist = getCosineContentDistriution(df2["PC2"])
PC3_wt_dist = getCosineContentDistriution(df2["PC3"])
plot(PC1_wt_dist)
plot(PC2_wt_dist)
plot(PC3_wt_dist)
xlabel("ts")
ylabel("Cosine content")
legend("PC1","PC2","PC3")
suptitle("Open wt cosine content of firsts 3 principal components")
show()"""