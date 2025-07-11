# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# By counting the steps between two adjacent pentagons (spiral centers),
# we find a path of 7 steps, a 60-degree turn, and then 2 steps.
# Given m > n, we have:
m = 7
n = 2

# Step 2: Determine the number of pentagonal faces (P).
# All Goldberg polyhedra have 12 pentagonal faces.
P = 12

# Step 3: Calculate the number of hexagonal faces (H).
# The formula is H = 10 * (T - 1), where T = m^2 + m*n + n^2.
T = m**2 + m*n + n**2
H = 10 * (T - 1)

# Step 4: Print the final answer in the format m,n,H,P.
print(f"{m},{n},{H},{P}")