import math

# This script calculates the exact value of 'c', the time of emergence
# of the giant connected component in a variant Erdos-Renyi random graph model.

# The analysis of the model shows that the critical point is reached when
# the average degree of the graph becomes 1. This leads to the following
# equation for the critical time c:
# c^2 / 3 = 1

# We will solve this equation programmatically.

# Step 1: Define the numbers in the final equation.
power = 2
divisor = 3
result = 1

print("The final equation for the critical time 'c' is:")
print(f"c^{power} / {divisor} = {result}")
print("-" * 30)

# Step 2: Solve the equation for c^2.
# c^2 = 3 * 1
c_squared = divisor * result
print("Solving for c^2:")
print(f"c^{power} = {divisor} * {result}")
print(f"c^{power} = {c_squared}")
print("-" * 30)

# Step 3: Solve for c by taking the square root.
# c = sqrt(3)
c_value = math.sqrt(c_squared)
print("Solving for c:")
print(f"c = sqrt({c_squared})")
print("-" * 30)

# Step 4: Output the final answer.
# The exact value is the square root of 3.
print(f"The exact value of c is sqrt({int(c_squared)})")
print(f"The numerical value of c is: {c_value}")
