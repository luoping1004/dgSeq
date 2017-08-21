import numpy as np
import scipy.optimize as op


def FindMin(Features, Label):
    nodeNum, featureNum = Features.shape
    
    def Gradient(theta, Features, Label):
        grad = np.zeros(featureNum)   
        for i in range(nodeNum):
            fw = 0
            fw = np.dot(theta,Features[i])
            for j in range(featureNum):
                grad[j] += Label[i]*Features[i,j] - (np.exp(fw)*Features[i,j])/(1+np.exp(fw))
        return -grad

    def ObjFunc(theta, Features, Label):
        #nodeNum, featureNum = Features.shape
        f = 0
        for i in range(nodeNum):
            fw = 0
            fw = np.dot(theta,Features[i])
            f += Label[i]*fw - np.log(1+np.exp(fw))
        return -f
    
    #nodeNum, featureNum = X.shape
    ini_theta = np.zeros(5)
    Result = op.minimize(fun = ObjFunc,
                         x0 = ini_theta,
                         args = (Features, Label),
                         method = 'Newton-CG',
                         jac = Gradient,
						 options={'maxiter':1000})
    return Result
