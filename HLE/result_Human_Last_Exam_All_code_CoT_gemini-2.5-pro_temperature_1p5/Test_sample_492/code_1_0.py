import math

# This script calculates the exact value of 'c', the time of emergence of the
# giant component in the described random graph model.

# The condition for the emergence of the giant component is that the
# average degree of the graph becomes 1. We will derive the average degree k(c)
# as a function of time c.

# Step 1: Find the expected number of vertices E[V(c)] and edges E[E(c)] at time c.
# - Vertices arrive at a stochastic rate n. This is a Poisson process.
#   E[V(c)] = n * c
#
# - Edges between any two existing vertices appear at rate 1/n. The expected rate
#   of new edge formation at time s is (1/n) * E[V(s) choose 2].
#   For a Poisson variable V(s) with mean n*s, E[V(s) choose 2] = (n*s)^2 / 2.
#   So, d/ds E[E(s)] = (1/n) * (n*s)^2 / 2 = n * s^2 / 2.
#
# - Integrating this rate from 0 to c gives the expected number of edges:
#   E[E(c)] = integral(n * s^2 / 2) ds from 0 to c = n * c^3 / 6.

# Step 2: Calculate the average degree k(c).
# k(c) is approximated by 2 * E[E(c)] / E[V(c)] for large n.
# k(c) = 2 * (n * c^3 / 6) / (n * c)
# k(c) = c^2 / 3

# Step 3: Set the average degree to 1 and solve for c.
# The equation is: k(c) = 1, which means c^2 / 3 = 1.
# This gives c^2 = 3.

# The numbers in the final equation k(c) = 1 are:
power = 2
denominator = 3
rhs = 1

# Solving c^2 = 3 for c (where c must be positive time)
c_squared_value = denominator * rhs
result = math.sqrt(c_squared_value)

print("The derivation of the average degree k(c) as a function of time c gives:")
print(f"k(c) = c^{power} / {denominator}")
print("\nThe giant component emerges when the average degree k(c) is 1.")
print("This gives the final equation to solve for the critical time c:")
print(f"c^{power} / {denominator} = {rhs}")
print(f"\nSolving this equation gives c^2 = {c_squared_value}, so the exact value is c = sqrt({c_squared_value}).")
print(f"The numerical value is c ≈ {result}")
