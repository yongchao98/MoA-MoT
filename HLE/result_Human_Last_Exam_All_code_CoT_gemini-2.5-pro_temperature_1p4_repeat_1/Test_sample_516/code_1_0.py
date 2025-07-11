import math

# The expected value E[T] is given by the formula: E[T] = 3 - 2 * e^(1/3).
# Here, we calculate this value.

# Define the constants in the equation
c1 = 3
c2 = -2
power = 1/3

# Calculate e^(1/3)
e_cubed_root = math.exp(power)

# Calculate the final result
result = c1 + c2 * e_cubed_root

# Print the equation with the calculated values
print(f"The calculation is based on the derived formula for the expected value E[T].")
print(f"E[T] = {c1} + ({c2}) * e^({power})")
print(f"e^({power}) is approximately {e_cubed_root}")
print(f"So, E[T] = {c1} - 2 * {e_cubed_root} = {result}")
