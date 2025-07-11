import math

# We are tasked with finding the time 'c' when a giant component emerges.
# This occurs when the average degree of the random graph equals 1.
# Let V(t) be the number of vertices and E(t) be the number of edges at time t.
# The condition for the giant component emergence is: 2 * E[E(c)] / E[V(c)] = 1

# From the problem description and mathematical derivation:
# 1. The expected number of vertices at time c is:
#    E[V(c)] = n * c
# 2. The expected number of edges at time c is:
#    E[E(c)] = (n * c^3) / 6

# Now, we set up the equation for the average degree being 1.
# 2 * ( (n * c^3) / 6 ) / ( n * c ) = 1

# Let's simplify this step by step:
# (n * c^3 / 3) / (n * c) = 1
# c^2 / 3 = 1

# The python code below solves this final simplified equation.

# Define the components of the simplified equation: c^2 / 3 = 1
lhs_denominator = 3
rhs_value = 1

# The equation is: c^2 / lhs_denominator = rhs_value
print("The equation for the critical time 'c' is derived from the average degree condition:")
print("2 * E[Number of Edges] / E[Number of Vertices] = 1")
print("")
print("Substituting the derived expectations gives:")
print("2 * (n * c^3 / 6) / (n * c) = 1")
print("")
print("This simplifies to the following equation:")
# We use f-string formatting to print the numbers in the final equation.
print(f"c^2 / {lhs_denominator} = {rhs_value}")
print("")

# Solving for c^2
c_squared = rhs_value * lhs_denominator
print(f"Solving for c^2, we get:")
print(f"c^2 = {c_squared}")
print("")

# Solving for c
c = math.sqrt(c_squared)
print(f"The exact value of c is the square root of {c_squared}:")
print(f"c = {c}")
