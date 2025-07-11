import math

# The problem reduces to finding the time `c` at which the average degree
# of the graph becomes 1.
# Based on the derivation, the average degree `d(t)` at time `t` is given by the formula:
# d(t) = t^2 / 3

# We need to solve the equation d(c) = 1.
# c^2 / 3 = 1

# Define the numbers in the final equation.
exponent = 2
divisor = 3
result = 1

print("The equation for the critical time c is:")
# Output each number in the final equation.
print(f"c^{exponent} / {divisor} = {result}")

# Solve for c^2
c_squared = result * divisor
print(f"\nFrom the equation, we find c^{exponent}:")
print(f"c^{exponent} = {c_squared}")

# Solve for c
# Since time must be positive, we take the positive square root.
c_value = math.sqrt(c_squared)
print("\nSolving for c (which must be positive):")
print(f"c = sqrt({c_squared})")

# Print the final exact value. Note: math.sqrt(3) returns a float.
# The exact value is the square root of 3.
print(f"\nThe exact value of c is: sqrt(3)")
print(f"The numerical value of c is approximately: {c_value}")