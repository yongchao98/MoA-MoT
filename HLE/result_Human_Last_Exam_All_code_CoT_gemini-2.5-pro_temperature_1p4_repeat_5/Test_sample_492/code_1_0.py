import math

# We are tasked with finding the time 'c' of the emergence of the giant
# connected component in a dynamic random graph.

# Our derivation shows that the average degree d(t) of the graph at time t
# follows the equation d(t) = t^2 / 3.
# The giant component emerges when the average degree becomes 1.
# This leads to the final equation for the critical time c.

# The equation has the form: c^2 / divisor = result
divisor = 3
result = 1

print(f"The equation for the critical time c is: c^2 / {divisor} = {result}")

# We solve for c:
# c^2 = divisor * result
c_squared = divisor * result
print(f"Solving for c^2, we get: c^2 = {c_squared}")

# c = sqrt(c_squared)
c = math.sqrt(c_squared)

print(f"The exact value of c is the square root of {c_squared}.")
print(f"The numerical value of c is: {c}")
