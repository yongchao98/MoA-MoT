import math

# Define the parameters of the problem.
# c: The propagation speed of the wave.
# T: The final time up to which convergence is required.
# M: The size of the overlap between subdomains.

# We will use some example values for demonstration.
c = 1.0
T = 5.0
M = 2.0

# The number of iterations N must satisfy N * (M / c) >= T.
# Rearranging for N gives N >= T * c / M.
# To find the smallest integer N, we use the ceiling function.

# Calculate the intermediate value before taking the ceiling.
value_before_ceil = T * c / M

# Calculate the final number of iterations.
num_iterations = math.ceil(value_before_ceil)

# Print the final result in the requested equation format.
# We cast the result to an integer for a cleaner output, as math.ceil() returns a float.
print(f"The required number of iterations N is calculated as:")
print(f"N = ceil(T * c / M)")
print(f"N = ceil({T} * {c} / {M})")
print(f"N = ceil({value_before_ceil})")
print(f"N = {int(num_iterations)}")