#################################################################
#-A script to perform autocorrelation analysis on a given array-#
#################################################################

#THIS SCRIPT IS WRONG#

def getProbability(value,distribution):
    """compute the probability of a vuale
    in a distribution"""
    n_apparence = distribution.count(value)
    return n_apparence/len(distribution)

def getExpectationValue(distribution):
    """compute the expectation value 
    of a distribution"""
    expectationValue = 0
    for value in distribution:
        expectationValue += value*getProbability(value,distribution)
    return expectationValue

def getVariance(distribution):
    """compute the variance 
    of a distibution"""
    variance = 0
    expectationValue = getExpectationValue(distribution)
    for value in distribution:
        variance += getProbability(value,distribution)*(value-expectationValue)**2
    return variance

def getAutocorrelationDistribution(distribution):
    """compute the autocorrelation ?function?
    of a distribution and rutern it as a dict"""
    from progressbar import ProgressBar
    distribution = list(distribution)
    expectationValue = getExpectationValue(distribution)
    autocorrelationDistribution = dict()
    variance = getVariance(distribution)
    progressbar = ProgressBar(max_value=len(distribution))
    for dt in range(len(distribution)):
        autoCorValueList = list()
        #print(dt)
        for counter in range(len(distribution)-dt):
            autoCorValue = (distribution[counter-1]-expectationValue)*(distribution[counter+dt-1]-expectationValue)
            autoCorValueList.append(autoCorValue)
        autocorrelationDistribution[dt] = getExpectationValue(autoCorValueList)/variance
        progressbar.update(dt)
    return autocorrelationDistribution


if __name__ == "__main__":
    from MDAnalysis import Universe
    from domains_angle import getAngle
    from pandas import Series
    from matplotlib.pyplot import plot,show,xlabel,ylabel,legend,suptitle


    #obj1 = Universe("MUTATED.pdb","final_whole_noPBC.xtc")
    obj2 = Universe("open_V_ps.pdb","open_450ns_protein_only.dcd")
    
    #angleResults1 = getAngle(obj1,roundOpt=True)
    angleResults2 = getAngle(obj2,sampling=2000,roundOpt=True)

    #angleResults1 = angleResults1[::10].reset_index(drop=True)



    #autoCorrelationResults1 = getAutocorrelationDistribution(list(angleResults1))
    autoCorrelationResults2 = getAutocorrelationDistribution(list(angleResults2))
    
    #toPlot1 = Series(autoCorrelationResults1)
    toPlot2 = Series(autoCorrelationResults2)

    plot(toPlot1,"r")
    plot(toPlot2,"b")
    xlabel("$[10^{-1}ns]$")
    ylabel("$C_f$")
    legend(("Mutated open conformation","Wild type open conformation"))
    suptitle("Autocorrelation of the angle between domains 1-2")

    show()
