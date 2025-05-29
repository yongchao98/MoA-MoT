import numpy as np
import math

def density(x, X):
    h = 0.25
    N = X.shape[0]
    accumulation = 0

    for x_i in X:
        accumulation += math.exp(-(x - x_i) ** 2 / (2 * (h ** 2))) / (h * math.sqrt(2 * math.pi))
    return accumulation / float(N)

# Test the hypothesis
X = np.array([0.0, 0.1, 0.2])
x_value = 15.0
density_value = density(x_value, X)
print(density_value)