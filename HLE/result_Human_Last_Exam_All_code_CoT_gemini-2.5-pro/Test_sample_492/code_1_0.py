import math

# This script calculates the constant 'c', the time of emergence of the giant
# connected component in the described variant of the Erdos-Renyi model.

# The derivation shows that the critical time 'c' is found by solving the
# equation k(c) = 1, where k(t) is the average degree of the graph at time t.
#
# The average degree k(t) is given by 2 * E(t) / V(t), where:
# - V(t) is the expected number of vertices ≈ n*t
# - E(t) is the expected number of edges ≈ n * t^3 / 6
#
# Substituting these into the condition k(c) = 1 gives:
# 2 * (n * c^3 / 6) / (n * c) = 1
# This simplifies to the final equation for c.

# The final equation to solve is: c^2 / 3 = 1
lhs_numerator_variable_power = 2
lhs_denominator = 3
rhs = 1

# Solving the equation step-by-step
c_squared = rhs * lhs_denominator
c = math.sqrt(c_squared)

print("The final equation for the critical time 'c' is derived from the condition that the average degree is 1.")
print(f"The equation is: c^{lhs_numerator_variable_power} / {lhs_denominator} = {rhs}")
print("\nSolving for c^2:")
print(f"c^2 = {rhs} * {lhs_denominator}")
print(f"c^2 = {c_squared}")
print("\nSolving for c:")
print(f"c = sqrt({c_squared})")
print(f"\nThe exact value of c is: {c}")
