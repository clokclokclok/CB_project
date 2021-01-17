from matplotlib.pyplot import plot,xlabel,ylabel,legend,suptitle,show


def xvgTOcsv(filename, write = False):
    """read .xvg GROMACS output files, if the parameter
    write is true it write also the .csv file"""
    from pandas import Series
    converted = list()
    with open(filename,"r") as toConvert:
        for line in toConvert:
            if not line.startswith("@") and not line.startswith("#"):
                line = line.strip(" ")
                line = line.strip("  ")
                line = line.strip("   ")
                line = line.strip("    ")
                line = line.replace("    ",",")
                line = line.replace("   ",",")
                line = line.replace("  ",",")
                line = line.replace(" ",",")
                converted.append(line)
    
    if write == True:
        new_fileName = filename.replace("xvg","csv")
        with open(new_fileName,"w") as toWrite:
            for line in converted:
                toWrite.write(line)
    
    resDict = dict()
    for line in converted:
        line = line.split(",")
        resDict[float(line[0])] = float(line[1])
    
    return Series(resDict)



if __name__ == "__main__" :

    file1 = "temperature.xvg"
    toPlot1 = xvgTOcsv(file1,write=True)
    runningMean1 = toPlot1.rolling(100).mean()
    plot(toPlot1,"b")
    plot(runningMean1,"r")
    xlabel("[ps]")
    ylabel("[K]")
    legend(("Temperatue","Averaged temperature(100 frames)"))
    suptitle("Temperature")
    show()

    file2 = "pressure.xvg"
    toPlot2 = xvgTOcsv(file2,write=True)
    runningMean2 = toPlot2.rolling(100).mean()
    plot(toPlot2,"b")
    plot(runningMean2,"r")
    xlabel("[ps]")
    ylabel("[bar]")
    legend(("Pressure","Averaged pressure(100 frames)"))
    suptitle("Pressure")
    show() 