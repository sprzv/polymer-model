import matplotlib.pyplot as plt
import numpy as np
from random import random
from scipy.optimize import curve_fit

#Chosen parameters

M = 4 #number of realisations of polymer chain
N = 400 #number of segments in chain

#Assign initial variables
histDataRG = np.zeros(M)

rEe2total = np.zeros(N + 1)
rG2total = np.zeros(N + 1)

histDataREe2 = np.zeros(M)
histDataRG2 = np.zeros(M)

histDataREe = np.zeros(M)

#figure for multiple single chain configurations
fig = plt.figure()

for k in range(1, M + 1): #k-th polymer configuration loop
    
    #reset all variables for single chain as new configuration loop starts
    rCmX = 0 #centre of mass x coordinate
    rCmY = 0 #centre of mass y coordinate
    
    xEnd = 0
    yEnd = 0 #current end of the chain in the x and y directions
    
    rX = np.zeros(N + 1)
    rY = np.zeros(N + 1)
    
    rX[0] = xEnd
    rY[0] = yEnd
    
    rCmXarray = np.zeros(N + 1)
    rCmYarray = np.zeros(N + 1)
    
    rXCom = np.zeros(N + 1)
    rYCom = np.zeros(N + 1)
    
    xArray = np.zeros(N + 1)
    yArray = np.zeros(N + 1)
    
    x1array = np.zeros(N + 1)
    y1array = np.zeros(N + 1)
    
    
    #Generate a polymer chain of N bonds (or N+1 sites)
    
    for i in range(1, N + 1):
        
        #add one bond at a time - generate coordinates for new bond point
        xB = random() - 0.5 #make sure new random numbers are generated for new chain generated
        yB = random() - 0.5
        
        #normalise length of displacement vector d so equal to 1 
        d = np.sqrt(xB**2 + yB**2) #subtract previous vector point though for each coord??
        x = xB / d
        y = yB / d
        
        #insert check of distance with other points here
        
        xArray[i] = x
        yArray[i] = y
        
        #new end of chain is at new bond point
        xEnd += x
        yEnd += y
        
        #note a new array focussed on where the mass is, is added
        rXCom[i] = xEnd - x/2
        rYCom[i] = yEnd - y/2
        
        #all the mass of a rod is placed 1/2 way along it
        rCmX += rXCom[i] #really the sum of all rCmX over i 
        rCmY += rYCom[i] #and same here
        
        rCmXarray[i] = rCmX
        rCmYarray[i] = rCmY
       
        #add new bond points into array of bond vectors
        rX[i] = xEnd
        rY[i] = yEnd
        
        #current end-to-end distance squared
        rEe2 = (rX[i] - rX[0])**2 + (rY[i] - rY[0])**2
        rEe = np.sqrt(rEe2)
        
        #add each rEe2 for each kth configuration in each ith array element slot - will be divided by M later to obtain average rEe2 as a function of N 
        #remember 0th element of array does not count - it's just 0 for N = 0
        rEe2total[i] = rEe2total[i] + rEe2
        
        
        #radius of gyration for i bonds in polymer chain
        rG2 = 0
        for j in range(1, i + 1):
            x1 = rXCom[j] - (rCmX / i) #centre of mass of current rod - centre of mass of whole chain
            y1 = rYCom[j] - (rCmY / i)
            
            x1array[j] = x1
            y1array[j] = y1
            
            rG2 += x1**2 + y1**2 #distance (squared) of site j + 1 (?) relative to centre of mass, sum over all rods i in the chain
        
        rG2 /= i  #should it be like this or divided by (i + 1)?? because need total number i = 0 to current i, so need to take into account the 0th element, which was not included in loop - he did divided by i
        rG = np.sqrt(rG2) #radius of gyration
        
        #add each rG2 for each kth configuration in each ith array element slot - will be divided by M later to obtain average rG2 as a function of N 
        #remember 0th element of array does not count - it's just 0 for N = 0
        rG2total[i] = rG2total[i] + rG2

        #print("i =", i)
        #print("rG =", rG)
        #print("rEe =", rEe)
        #print("rEe / rG =", rEe / rG)
        
    #end of bond addition loop
    
    plt.plot(x1array, y1array, '-', label = 'k = {}'.format(k))
    plt.plot(x1array[0], y1array[0], '^')
    plt.plot(x1array[N], y1array[N], 'v')
    
    histDataREe2[k - 1] = rEe2 #takes last rEe2 calculated, the one for i = N
    histDataRG2[k - 1] = rG2
    
    histDataREe[k - 1] = rEe #takes last rEe2 calculated, the one for i = N
    histDataRG[k - 1] = rG
    
    #print("rEe2total for M = {}".format(k), "is", rEe2total)
    #print("rG2total for M = {}".format(k), "is", rG2total)


#end of polymer configuration loop

#outside all loops - at end of program
rEeAverage = np.zeros(N + 1)
rGAverage = np.zeros(N + 1)
rGrEeRatio = np.zeros(N + 1)

for l in range(1, N + 1):
    rEe2total[l] = rEe2total[l] / M
    rG2total[l] = rG2total[l] / M
    rEeAverage[l] = np.sqrt(rEe2total[l])
    rGAverage[l] = np.sqrt(rG2total[l])
    rGrEeRatio[l] = rGAverage[l] / rEeAverage[l]

#print("rEeAverage =", rEeAverage)
#print("rGAverage =", rGAverage)

#plot example displacement for multiple single chain configurations

plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.legend()


#plotting rEe2 and rG2 as a function of N - plot in same graph ?? what values of M for each N are needed for reliable statistics??
arrayN = np.arange(0, N + 1)
xPlot = np.linspace(0, N + 1, 1000)
yPlot = xPlot

fig1 = plt.figure()
plt.plot(arrayN, rEe2total, label = 'N = {}, M = {}'.format(N, M))
plt.plot(xPlot, yPlot, 'r', label = 'y = x')
plt.xlabel("N")
plt.ylabel("$R_{ee}^2$")
plt.legend()

fig2 = plt.figure()
plt.plot(arrayN, rG2total, label = 'N = {}, M = {}'.format(N, M))
plt.plot(xPlot, yPlot / 6, 'r',label = 'y = 1/6 x')
plt.xlabel("N")
plt.ylabel("$R_{g}^2$")
#plt.xscale('log')
#plt.yscale('log')
plt.legend()

##log log plots
#fig1 = plt.figure()
#plt.plot(arrayN, rEeAverage, label = 'N = {}, M = {}'.format(N, M))
#plt.xlabel("N")
#plt.ylabel("Average $R_{ee}$")
#plt.xscale('log')
#plt.yscale('log')
#plt.legend()
#
#fig2 = plt.figure()
#plt.plot(arrayN, rGAverage, label = 'N = {}, M = {}'.format(N, M))
#plt.xlabel("N")
#plt.ylabel("Average $R_{g}$")
#plt.xscale('log')
#plt.yscale('log')
#plt.legend()


#plot ratio of average Rg and Ree as function of N
fig3 = plt.figure()
plt.plot(arrayN, rGrEeRatio, label = 'M = {}'.format(N, M))
plt.xlabel("N")
plt.ylabel("Average $R_{g}$ / Average $R_{ee}$")
#plt.xscale('log')
#plt.yscale('log')
plt.legend()


#use fit functions to find exact power for N?
def func(x, *theta):
    nu, alpha, beta = theta
    polyfunc = beta + alpha* x**(nu)
    return polyfunc


#sig = np.zeros(N + 1)
###random sigma should actually calculate it to obtain accurate nu result
#for i in range(0, N+1):
#    sig[i] = np.sqrt(i)
#    print(sig[i])
#    
### set default parameter values and do the fit
#p0 = np.array([1.0, 1.0, 0.0])
#thetaHat, cov = curve_fit(func, arrayN, rGAverage, p0, sig, absolute_sigma=True)
#print(thetaHat[0])

#histograms for fixed N
#fig3 = plt.figure()
#plt.hist(histDataREe2, bins = 50)
#plt.xlabel("$R_{ee}^2$")
#plt.ylabel("Counts")

#fig4 = plt.figure()
#plt.hist(histDataRG2, bins = 50)
#plt.xlabel("$R_{g}^2$")
#plt.ylabel("Counts")

fig3 = plt.figure()
plt.hist(histDataREe, bins = 40)
plt.xlabel("$R_{ee}$")
plt.ylabel("Counts")

fig4 = plt.figure()
plt.hist(histDataRG, bins = 40)
plt.xlabel("$R_{g}$")
plt.ylabel("Counts")

plt.show()