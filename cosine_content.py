from math import cos,pi
from numba import jit

@jit
def getCosineContent(distribution,steps=None,cosConst = 1):
    if steps == None:
        steps = len(distribution)
    integrationOfNumerator = 0
    integrationOfDenominator = 0
    for step in range(steps):
        integrationOfNumerator += cos(pi*cosConst*step)*distribution[step]
        integrationOfDenominator += distribution[step]**2
    cosineContent = integrationOfNumerator**2/integrationOfDenominator*(2/steps)
    return cosineContent

def getCosineContentDistriution(distribution,resolution=1,cosConst = 1):
    from progressbar import ProgressBar
    from pandas import Series
    bar = ProgressBar(max_value=len(distribution))
    cosResDist = {}
    for ts in range(resolution,len(distribution),resolution):
        cosResDist[ts] = getCosineContent(distribution,steps=ts,cosConst=cosConst)
        bar.update(ts)
    return Series(cosResDist)

 
if __name__ == "__main__" :
    from pandas import DataFrame,Series,read_csv
    from matplotlib.pyplot import plot,show,xlabel,ylabel,suptitle,legend
    from os import mkdir,listdir


    #importing PCA's saved results
    df2 = read_csv("Mutated_pca.csv")
    df1 = read_csv("Mutated_closed_pca.csv")


    #lowering resolutions of df1 series
    pc1 = df1["PC1"]
    pc2 = df1["PC2"]
    pc3 = df1["PC3"]
    low_pc1 = pc1[::10]
    low_pc1 = low_pc1.reset_index(drop = True)
    low_pc2 = pc2[::10]
    low_pc2 = low_pc2.reset_index(drop = True)
    low_pc3 = pc3[::10]
    low_pc3 = low_pc3.reset_index(drop = True)

    #lowering resolutions of df2 series
    pc1_2 = df2["PC1"]
    pc2_2 = df2["PC2"]
    pc3_2 = df2["PC3"]
    low2_pc1 = pc1_2[::10]
    low2_pc1 = low2_pc1.reset_index(drop = True)
    low2_pc2 = pc2_2[::10]
    low2_pc2 = low2_pc2.reset_index(drop = True)
    low2_pc3 = pc3_2[::10]
    low2_pc3 = low2_pc3.reset_index(drop = True)
    

    #computing and saving cosine content distribution per ts for PCs
    if "cosine_content_results" not in listdir():
        mkdir("cosine_content_results")
    PC1_mutated_dist = getCosineContentDistriution(low_pc1,cosConst=1/2)
    PC1_mutated_dist.to_csv("./cosine_content_results/pc1_mutated_closed.csv")
    PC2_mutated_dist = getCosineContentDistriution(low_pc2)
    PC2_mutated_dist.to_csv("./cosine_content_results/pc2_mutated_closed.csv")
    PC3_mutated_dist = getCosineContentDistriution(low_pc3,cosConst=3/2)
    PC3_mutated_dist.to_csv("./cosine_content_results/pc3_mutated_closed.csv")
    plot(PC1_mutated_dist)
    plot(PC2_mutated_dist)
    plot(PC3_mutated_dist)
    xlabel("frame(reduction factor of 10)")
    ylabel("Cosine content")
    legend(("PC1","PC2","PC3"))
    suptitle("Closed mutated cosine content of firsts 3 principal components")
    show()

    PC1_wt_dist = getCosineContentDistriution(low2_pc1,cosConst=1/2)
    PC1_wt_dist.to_csv("./cosine_content_results/pc1_mutated.csv")
    PC2_wt_dist = getCosineContentDistriution(low2_pc2)
    PC2_wt_dist.to_csv("./cosine_content_results/pc2_mutated.csv")
    PC3_wt_dist = getCosineContentDistriution(low2_pc3,cosConst=3/2)
    PC3_wt_dist.to_csv("./cosine_content_results/pc3_mutated.csv")
    plot(PC1_wt_dist)
    plot(PC2_wt_dist)
    plot(PC3_wt_dist)
    xlabel("frame(reduction factor of 10)")
    ylabel("Cosine content")
    legend(("PC1","PC2","PC3"))
    suptitle("Open mutated cosine content of firsts 3 principal components")
    show()