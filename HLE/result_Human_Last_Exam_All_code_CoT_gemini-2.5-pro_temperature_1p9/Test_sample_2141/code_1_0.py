import math

# This script calculates the maximum value of the ratio D_n(r*)/D_n^c(r*)
# As derived in the explanation, the maximum occurs at n=1.
# The maximum value is given by the expression: 4 * sqrt(3) * pi * e^(-3)

# Define the components of the final equation
val_4 = 4
sqrt_3 = math.sqrt(3)
pi = math.pi
exp_neg_3 = math.exp(-3)

# Calculate the final result
result = val_4 * sqrt_3 * pi * exp_neg_3

# Print the components and the final answer
print("The final expression for the maximum ratio is: 4 * sqrt(3) * pi * exp(-3)")
print(f"The value of each component is:")
print(f"4 = {val_4}")
print(f"sqrt(3) = {sqrt_3}")
print(f"pi = {pi}")
print(f"exp(-3) = {exp_neg_3}")
print("\nThe final numerical value is:")
print(result)
