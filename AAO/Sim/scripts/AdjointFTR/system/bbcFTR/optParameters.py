import numpy as np

def getParamArray(paramObj):
    arr = []
    for i,elem in enumerate(paramObj):
        for elemName in elem:
            arr.append(paramObj[i][elemName])
    return np.array(arr)

def getParamObj(paramArray):
    params = []
    numParam = 3
    numMags = int(len(paramArray) / numParam)
    for i in range(numMags):
        tmp = {}
        tmp['zstart'] = paramArray[i*numParam]
        tmp['rotation'] = paramArray[i*numParam+1]
        tmp['dbdx'] = paramArray[i*numParam+2]
        params.append(tmp)
    return np.array(params)   

if False:
    paramarray = np.array([
    1.157924056225655,
    -0.977797295545406,
    0.999016701746482,
    1.015198307362092,
    1.021911426027822,
    1.040902706986164,
    1.009810959559711,
    1.044424737861931,
    1.50542061709607,
    1.030422542262641,
    1.021398090036066 ,   
    ])

    #paramarray = np.ones(11)

    quad1 = {}
    quad1['zstart'] = paramarray[2]
    quad1['rotation'] = paramarray[8]
    quad1['dbdx'] = paramarray[3]

    quad2 = {}
    quad2['zstart'] = paramarray[4]
    quad2['rotation'] = paramarray[9]
    quad2['dbdx'] = paramarray[5]

    quad3 = {}
    quad3['zstart'] = paramarray[6]
    quad3['rotation'] = paramarray[10]
    quad3['dbdx'] = paramarray[7]

    sol = {}
    sol['zstart'] = paramarray[0]
    sol['rotation'] = 1.0
    sol['dbdx'] = paramarray[1]

if True:

    paramarray = np.array([ 
        0.97349958,  
        1.17730812,  
        1.22577242,  
        1.6338717,   
        1.13054699,  
        1.31214369,
        1.51125034,  
        1.09008688,  
        1.01156225,  
        0.8,         
        0.99995458, 
        -0.94121212
    ])    

    quad1 = {}
    quad1['zstart'] = paramarray[0]
    quad1['rotation'] = paramarray[1]
    quad1['dbdx'] = paramarray[2]

    quad2 = {}
    quad2['zstart'] = paramarray[3]
    quad2['rotation'] = paramarray[4]
    quad2['dbdx'] = paramarray[5]

    quad3 = {}
    quad3['zstart'] = paramarray[6]
    quad3['rotation'] = paramarray[7]
    quad3['dbdx'] = paramarray[8]

    sol = {}
    sol['zstart'] = paramarray[9]
    sol['rotation'] = 1.0
    sol['dbdx'] = paramarray[11]    

params = np.array([ quad1, quad2, quad3, sol])

# Opt parameters
# quad1 = {}
# quad1['zstart'] = 0.999016701746482
# quad1['rotation'] = 1.030542061709607
# quad1['dbdx'] = 1.015198307362092

# quad2 = {}
# quad2['zstart'] = 1.021911426027822
# quad2['rotation'] = 1.030422542262641
# quad2['dbdx'] = 1.040902706986164

# quad3 = {}
# quad3['zstart'] = 1.009810959559711
# quad3['rotation'] = 1.021398090036066
# quad3['dbdx'] = 1.044424737861931

# sol = {}
# sol['zstart'] = 1.157924056225655
# sol['rotation'] = 1.0
# sol['dbdx'] = -0.977797295545406

params = np.array([ quad1, quad2, quad3, sol])