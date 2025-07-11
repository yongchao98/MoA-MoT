# Parameters identified from visual inspection of the Goldberg polyhedron.
# m represents the number of steps in a straight line between adjacent pentagons.
m = 4
# n represents the number of steps after a 60-degree turn.
n = 2

# For a Goldberg polyhedron derived from an icosahedron, the number of pentagonal faces is always 12.
P = 12

# The number of hexagonal faces (H) is calculated using the formula:
# H = 10 * (T - 1), where T = m^2 + m*n + n^2.
T = m**2 + m * n + n**2
H = 10 * (T - 1)

# The problem asks for the answer in the format m,n,H,P without spaces.
# The following print statement assembles the final string.
# Each number in the final answer is explicitly represented.
print(f"{m},{n},{H},{P}")