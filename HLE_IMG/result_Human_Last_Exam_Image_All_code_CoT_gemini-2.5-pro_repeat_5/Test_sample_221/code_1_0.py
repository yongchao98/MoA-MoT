# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# From counting the faces between two adjacent pentagons, we found m + n = 6.
# Given the hint m > n and the chiral nature (n > 0), we determined (m, n) = (4, 2).
m = 4
n = 2

# Step 2: Calculate the number of pentagonal and hexagonal faces.
# For any Goldberg polyhedron based on an icosahedron, the number of pentagonal faces is 12.
P = 12

# The number of hexagonal faces H is calculated using the triangulation number T.
# T = m^2 + m*n + n^2
T = m**2 + m * n + n**2

# H = 10 * (T - 1)
H = 10 * (T - 1)

# Step 3: Print the final answer in the format m,n,H,P.
# The problem asks to output each number in the final equation. We will print the final values.
print(f"{m},{n},{H},{P}")
