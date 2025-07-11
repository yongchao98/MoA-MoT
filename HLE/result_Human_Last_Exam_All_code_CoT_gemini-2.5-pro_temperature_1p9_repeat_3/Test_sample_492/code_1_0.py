import sympy

# This script calculates the exact time 'c' of the emergence of the
# giant connected component in the specified random graph model.

# Step 1: State the derived formulas for the expected number of vertices and edges
# at the critical time 'c'. In the large n limit:
# E[N(c)] = n * c
# E[M(c)] = n * c^3 / 6
print("Let 'c' be the time of the emergence of the giant component.")
print("The expected number of vertices at time c is E[N(c)] = n * c.")
print("The expected number of edges at time c is E[M(c)] = n * c^3 / 6.")
print("-" * 30)

# Step 2: Apply the critical condition for the giant component.
# This occurs when the average degree of the graph is equal to 1.
# Average Degree = 2 * E[M(c)] / E[N(c)]
print("The giant component emerges when the average degree equals 1:")
print("2 * E[M(c)] / E[N(c)] = 1")
print("2 * (n * c^3 / 6) / (n * c) = 1")
print("-" * 30)


# Step 3: Simplify and present the final equation for 'c'.
# The parameter 'n' cancels out in the simplification process.
# => 2 * c^3 / (6 * c) = 1
# => c^2 / 3 = 1
print("The final equation for c is:")
# For the equation c^2 / 3 = 1, the numbers are 2 (as a power), 3, and 1.
print("c^2 / 3 = 1")
print("-" * 30)

# Step 4: Solve the equation to find the exact value of c.
# Define 'c' as a positive symbolic variable to get the correct root.
c = sympy.Symbol('c', positive=True)

# Define the equation c**2 / 3 = 1
equation = sympy.Eq(c**2 / 3, 1)

# Solve the equation for c
solution = sympy.solve(equation, c)

# The result is a list, so we extract the single element.
final_c_value = solution[0]

print("Solving the equation for c:")
print("c^2 = 3")
print(f"The exact value of c is {final_c_value}")
