#################################################################
#-A script to perform autocorrelation analysis on a given array-#
#################################################################
from numba import jit

@jit(nopython=True)
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

@jit(nopython=True)
def getAutocorrelationDistributionINTERNAL(distribution,percentage,expectationValue,variance):
    autocorrelationDistribution = dict()
    for dt in range(len(distribution)):
        print("Computing autocorrelatione function:{}".format(dt/percentage))
        autoCorValueList = list()
        for counter in range(len(distribution)-dt):
            autoCorValue = (distribution[counter-1]-expectationValue)*(distribution[counter+dt-1]-expectationValue)
            autoCorValueList.append(autoCorValue)
        autocorrelationDistribution[dt] = getExpectationValue(autoCorValueList)/variance
    return autocorrelationDistribution
    
def getAutocorrelationDistribution(distribution):
    """compute the autocorrelation ?function?
    of a distribution and rutern it as a dict"""
    distribution = list(distribution)
    expectationValue = getExpectationValue(distribution)
    percentage = len(distribution)/100
    variance = getVariance(distribution)
    autoCorrelationResults = getAutocorrelationDistributionINTERNAL(distribution,percentage,expectationValue,variance)
    return autoCorrelationResults

if __name__ == "__main__":
    from MDAnalysis import Universe
    from domains_angle import getAngle
    from pandas import Series
    from matplotlib.pyplot import plot,show

    obj = Universe("open_V_ps.pdb","open_450ns_protein_only.dcd")
    angleResults = getAngle(obj,sampling=2000)
    autoCorrelationResults = getAutocorrelationDistribution(list(angleResults))

    toPlot = Series(autoCorrelationResults)

    plot(toPlot)
    show()