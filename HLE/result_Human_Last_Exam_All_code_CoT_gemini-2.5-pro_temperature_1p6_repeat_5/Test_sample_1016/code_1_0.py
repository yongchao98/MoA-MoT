import math

# --- Problem Parameters ---

# Final time of the simulation
T = 50.0

# Propagation speed of the wave
c = 2.0

# Left boundary of the right subdomain
a = 10.0

# Right boundary of the left subdomain
b = 12.5

# --- Calculation ---

# Calculate the size of the overlap region
M = b - a

# The formula for the number of iterations N is floor(T * c / M).
# This is derived from the time it takes for information to propagate
# back and forth across the subdomains until the entire spacetime
# domain [0, L] x [0, T] is correctly resolved.

# We need K > T*c/M total updates.
# The minimum integer K is floor(T*c/M) + 1.
# The counter N is K-1.
# So, N = floor(T*c/M).
num_iterations = math.floor(T * c / M)


# --- Output ---

# The final output should show the equation with all the numbers.
print(f"Given T = {T}, c = {c}, and an overlap M = b - a = {b} - {a} = {M}:")
print("The required number of iterations (N) is calculated as:")
print(f"N = floor(T * c / M)")
print(f"N = floor({T} * {c} / {M})")
# Intermediate step to show the value inside floor
value_inside_floor = T * c / M
print(f"N = floor({value_inside_floor})")
print(f"N = {num_iterations}")

# The final answer in the requested format will be the last calculated value.
# For T=50, c=2, M=2.5, N = floor(50 * 2 / 2.5) = floor(100 / 2.5) = floor(40) = 40.