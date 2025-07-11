import numpy as np

# Let's assume a value for N for demonstration purposes
N = 10
# The derived factor for the upper bound
factor = 2

# The upper bound is factor * sqrt(N)
upper_bound_expression = f"{factor} * sqrt(N)"

# Let's print the parts of the final equation
# The question asks to output each number in the final equation.
# Here, the only number derived is the factor '2'.
print("The upper bound for ||B Q_{0,M}||_oo is given by K * sqrt(N)")
print(f"The derived value for the constant factor K is:")
print(factor)
print("\nFor a sample graph with N = 10 nodes, the equation for the bound would be:")
print(f"Bound = {factor} * sqrt({N}) = {factor * np.sqrt(N):.2f}")