# The parameter d is a positive even integer.
# The number of edges to add in the worst case for a given d is K(d) = 3*d/2 + 1.
# The problem asks for the minimal possible value of this quantity.
# Since K(d) is an increasing function of d, the minimum occurs at the smallest possible value for d, which is d=2.

d = 2
minimal_number_of_edges = 3 * d / 2 + 1

# Output the equation with the numbers plugged in.
print(f"The minimal number of edges is given by the formula 3*d/2 + 1 for the smallest possible d.")
print(f"The smallest even degree d for a 2-edge-connected graph is 2.")
print(f"Plugging in d = {d}:")
print(f"3 * {d} / 2 + 1 = {int(minimal_number_of_edges)}")