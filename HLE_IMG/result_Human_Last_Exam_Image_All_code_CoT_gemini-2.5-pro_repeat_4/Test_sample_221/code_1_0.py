# Step 1: Define the parameters found from visual inspection.
# By counting the hexagons between two adjacent pentagons, we found a path of 4 hexagons.
# For a straight path (n=0), the number of hexagons is m-1.
# So, m - 1 = 4, which means m = 5.
m = 5
n = 0

# Step 2: The number of pentagonal faces in a Goldberg polyhedron is always 12.
P = 12

# Step 3: Calculate the number of hexagonal faces using the formula.
# H = 10 * (m^2 + m*n + n^2 - 1)
H = 10 * (m**2 + m*n + n**2 - 1)

# Step 4: Print the final answer in the required format "m,n,H,P".
# The output should not have spaces.
print(f"{m},{n},{H},{P}")