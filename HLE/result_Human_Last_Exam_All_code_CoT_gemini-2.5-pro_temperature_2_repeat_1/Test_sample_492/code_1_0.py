import math

# This script calculates the critical time 'c' for the emergence of a giant
# component in the described random graph model.

# The derivation shows that the average degree of the graph at time 't' is d(t) = t^2 / 3.
# The giant component emerges when the average degree d(c) = 1.
# This gives us the equation: c^2 / 3 = 1.

# We define the numbers from the final equation.
# In the equation c^2 / divisor = result, we have:
divisor = 3
result = 1

# We print the equation with its numerical components.
print(f"The final equation for the critical time 'c' is: c^2 / {divisor} = {result}")

# To solve for c, we first solve for c^2.
# c^2 = result * divisor
c_squared = result * divisor
print(f"Rearranging the equation gives: c^2 = {result} * {divisor}")
print(f"This simplifies to: c^2 = {c_squared}")


# Then, we take the square root to find c.
c = math.sqrt(c_squared)

print(f"Finally, solving for c gives: c = sqrt({c_squared})")
print(f"\nThe exact value of c is the square root of 3.")
print(f"The numerical value of c is approximately: {c}")