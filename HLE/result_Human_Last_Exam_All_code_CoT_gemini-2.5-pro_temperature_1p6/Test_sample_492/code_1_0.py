import math

# This script calculates the exact value of 'c', the time of emergence of the
# giant component in the described graph model.

# The condition for emergence is when the average degree k(c) = 1.
# The average degree k(t) as a function of time t is given by k(t) = t^2 / 3.
# This leads to the final equation for c: c^2 / 3 = 1.

# We define the numbers from the final equation: c^exponent / denominator = result
exponent = 2
denominator = 3
result = 1

# Print the final equation with its components
print("The final equation for the critical time 'c' is derived from the condition that the average degree is 1:")
print(f"c^{exponent} / {denominator} = {result}")
print("-" * 30)

# Solve the equation for c
# c^2 = 3
# c = sqrt(3)
c_squared = result * denominator
print("Solving for c:")
print(f"c^{exponent} = {c_squared}")

c = math.sqrt(c_squared)
print(f"c = sqrt({c_squared})")
print("-" * 30)

# Print the final answer
print(f"The exact value of c is sqrt({c_squared}).")
print(f"The numerical value of c is: {c}")