import numpy as np

def fx(inX):
    return inX**3 - 6*inX**2 + 11*inX - 6

def deriv_fx(inX):
    return 3*inX**2 - 12*inX + 11

def GradientAscent(startX, endX, alpha=0.02):
    outX = []
    xOld = startX - 1
    xNew = startX
    epsilon = 0.000001
    while xNew - xOld > epsilon and xNew < endX:
        xOld = xNew
        outX.append(xOld)
        xNew = xOld + alpha * deriv_fx(xOld)
    outX = np.array(outX)
    outY = fx(outX)
    return outX, outY

# Further adjusted inputs
startX = 0.993
endX = 1.5
alpha = 0.02

outX, outY = GradientAscent(startX, endX, alpha)
print({"x_path": outX.tolist(), "y_path": outY.tolist()})