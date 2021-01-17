from matplotlib.pyplot import plot,xlabel,ylabel,legend,suptitle,show
from os import chdir
from pandas import read_csv


chdir("cosine_content_results")

#plotting mutated closed
toPlot1 = read_csv("pc1_mutated_closed.csv",index_col=0,squeeze=True)
toPlot2 = read_csv("pc2_mutated_closed.csv",index_col=0,squeeze=True)
toPlot3 = read_csv("pc3_mutated_closed.csv",index_col=0,squeeze=True)
plot(toPlot1)
plot(toPlot2)
plot(toPlot3)
ylabel("Cosine content")
xlabel("[$10^1$ns]")
legend(("PC1","PC2","PC2"))
suptitle("Closed mutated cosine content of firsts 3 principal components")
show()

#plotting mutated open
toPlot1 = read_csv("pc1_mutated.csv",index_col=0,squeeze=True)
toPlot2 = read_csv("pc2_mutated.csv",index_col=0,squeeze=True)
toPlot3 = read_csv("pc3_mutated.csv",index_col=0,squeeze=True)
plot(toPlot1)
plot(toPlot2)
plot(toPlot3)
ylabel("Cosine content")
xlabel("[$10^1$ns]")
legend(("PC1","PC2","PC2"))
suptitle("Open mutated cosine content of firsts 3 principal components")
show()

#plotting wt open
toPlot1 = read_csv("pc1_wt.csv",index_col=0,squeeze=True)
toPlot2 = read_csv("pc2_wt.csv",index_col=0,squeeze=True)
toPlot3 = read_csv("pc3_wt.csv",index_col=0,squeeze=True)
plot(toPlot1)
plot(toPlot2)
plot(toPlot3)
ylabel("Cosine content")
xlabel("[$10^1$ns]")
legend(("PC1","PC2","PC2"))
suptitle("Open wild type cosine content of firsts 3 principal components")
show()