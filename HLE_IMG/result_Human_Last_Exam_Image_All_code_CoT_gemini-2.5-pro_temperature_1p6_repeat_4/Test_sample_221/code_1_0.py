# Step 1: Define the parameters m and n based on visual inspection.
# By tracing the path between two adjacent pentagons along the hexagonal grid,
# we can count the steps. We move m steps, turn 60 degrees, and move n steps.
# Visual analysis of the image reveals a path of 4 steps and 2 steps.
# As per the hint, m > n.
m = 4
n = 2

# Step 2: Define the number of pentagonal faces (P).
# All Goldberg polyhedra have 12 pentagonal faces.
P = 12

# Step 3: Calculate the number of hexagonal faces (H).
# The formula for the number of hexagonal faces is H = 10 * (m^2 + m*n + n^2 - 1).
# We calculate the value of T = m^2 + m*n + n^2 first.
T = m**2 + m*n + n**2
H = 10 * (T - 1)

# Step 4: Print the final answer in the specified format "m,n,H,P".
print(f"{m},{n},{H},{P}")