import numpy as np
from scipy.integrate import quad

# The final temperature T is given by the equation:
# T = -2 + e^4 / (e + I)
# where I is the definite integral of exp(x^2) from 1 to 2.

# 1. Define the function to be integrated
def integrand(x):
    return np.exp(x**2)

# 2. Calculate the definite integral I = integral from 1 to 2 of exp(x^2) dx
# quad returns a tuple (integral_value, estimated_error)
I, error = quad(integrand, 1, 2)

# 3. Get the values for e and e^4
e = np.e
e4 = np.exp(4)

# 4. Calculate the final temperature T at the target time
T_final = -2 + e4 / (e + I)

# 5. Print the equation with the computed numerical values
print(f"The equation for the final temperature is:")
print(f"T = -2 + exp(4) / (exp(1) + integral from 1 to 2 of exp(x^2) dx)")
print(f"Substituting the numerical values:")
# We use round() for a cleaner display of the numbers
print(f"T = -2 + {round(e4, 5)} / ({round(e, 5)} + {round(I, 5)})")

# 6. Print the final result
print("\nThe final temperature at the specified time is:")
print(T_final)