# (a) Does the property of piecewise polynomiality of Z imply continuity?
# As explained in the thinking steps, the answer is No. A function can be
# piecewise polynomial without being continuous (e.g., a step function).

# (b) For the specific case g=0, n_+=3, n_-=1, determine the degree.
# We will calculate this using the formula for the degree of the polynomial Z(g,n).

# Parameters from the problem
g = 0
n_plus = 3
n_minus = 1

# The total number of boundaries, n, is the sum of n_+ and n_-.
n = n_plus + n_minus

# The degree of the polynomial Z(g, n) is given by the formula 3g - 3 + n.
degree = 3 * g - 3 + n

print("Part (a):")
print("No, the property of being piecewise polynomial does not by itself imply continuity.")
print("-" * 20)
print("Part (b):")
print(f"To find the degree of the polynomial Z for g={g}, n_+={n_plus}, n_-={n_minus}, we use the formula: Degree = 3g - 3 + n.")
print(f"First, we calculate the total number of boundaries, n:")
print(f"n = n_+ + n_- = {n_plus} + {n_minus} = {n}")
print("\nNow, we substitute g and n into the degree formula:")
print(f"Degree = 3 * g - 3 + n")
print(f"Degree = 3 * {g} - 3 + {n}")
print(f"Degree = {3*g} - 3 + {n}")
print(f"Degree = {3*g - 3} + {n}")
print(f"Degree = {degree}")
