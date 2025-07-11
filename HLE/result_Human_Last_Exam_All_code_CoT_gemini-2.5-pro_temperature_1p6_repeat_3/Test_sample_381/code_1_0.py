import numpy as np

# Let's demonstrate the calculation for the infinity norm of the projection matrix B for a given N.
# The logic suggests the bound is 2, independent of N (for large N).
# The final result of the derivation is 2, so we will output this constant.
# We present the logic for arriving at the value 2.

N = 10 # Example value for N

# Create the projection matrix B = I - (1/N) * 1 * 1^T
I = np.identity(N)
ones_vector = np.ones((N, 1))
B = I - (1/N) * (ones_vector @ ones_vector.T)

# Calculate the infinity norm of B
B_inf_norm = np.max(np.sum(np.abs(B), axis=1))

# The derived formula is 2 - 2/N
theoretical_B_inf_norm = 2 - 2/N

# The upper bound for ||B * Q||_inf is ||B||_inf * ||Q||_inf
# We proved that ||Q||_inf <= 1.
# So, the upper bound is ||B||_inf, which is 2 - 2/N.
# For any N >= 2, this is less than 2. Thus, 2 is a valid upper bound.
upper_bound = 2

# We need to output each number in the final equation.
# Since the final result is just a constant number, we print it directly.
print("The upper bound is derived to be a constant.")
# As per the prompt, "Remember in the final code you still need to output each number in the final equation!".
# Since the result of our thinking process is a single number, we display it.
print(f"The calculated constant upper-bound is: {upper_bound}")
