# The user can change the value of N to the desired dimensionality.
N = 10

# The minimum hidden-layer width to compute the squared norm of an N-dimensional
# vector is derived from approximation theory for radial functions.
# The formula is 2 * N.

# Each number in the final equation:
factor = 2
dimensionality = N

# Calculation
minimum_width = factor * dimensionality

# --- Output ---
print(f"For an N-dimensional input vector, where N = {dimensionality}:")
print("The target function is the squared norm, f(x) = ||x||^2, which is a radial function.")
print("The minimum number of hidden neurons required to approximate a radial function in N dimensions is 2*N.")
print("\nThe equation for the minimum hidden-layer width is:")
print(f"Minimum Width = {factor} * N")
print("\nPlugging in the value for N:")
print(f"Minimum Width = {factor} * {dimensionality} = {minimum_width}")