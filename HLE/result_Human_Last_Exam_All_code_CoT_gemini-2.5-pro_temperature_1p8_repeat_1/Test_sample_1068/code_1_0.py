import numpy as np

# This script calculates the value of the function l(p) at p=14.
# After simplification, the function is found to be l(p) = 2 * p * log(2).
# For p = 14, this becomes l(14) = 2 * 14 * log(2) = 28 * log(2).

# Define the numbers in the final simplified equation l(14) = coeff * log(arg)
p = 14
coeff = 2 * p
arg = 2

# Calculate the result
result = coeff * np.log(arg)

# Print the final equation with its components
print(f"The simplified symbolic expression for l({p}) is {coeff}*log({arg}).")

# Print the final numerical value
print(f"The numerical value of l({p}) is: {result}")