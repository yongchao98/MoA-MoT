# Parameters of the Goldberg Polyhedron identified from the image.
# m represents the number of steps in a straight line.
# n represents the number of steps after a 60-degree turn.
m = 5
n = 3

# The number of pentagonal faces (P) in a Goldberg polyhedron is always 12.
P = 12

# Calculate the triangulation number (T).
# The equation is T = m^2 + m*n + n^2
T = m**2 + m * n + n**2

# Calculate the number of hexagonal faces (H).
# The equation is H = 10 * (T - 1)
H = 10 * (T - 1)

# Print the final result in the format m,n,H,P as requested.
# We also show the calculation as per the instructions.
print(f"The parameters are m={m} and n={n}.")
print(f"The number of pentagonal faces P is {P}.")
print(f"The calculation for hexagonal faces H is:")
print(f"T = {m}^2 + ({m} * {n}) + {n}^2 = {m**2} + {m*n} + {n**2} = {T}")
print(f"H = 10 * (T - 1) = 10 * ({T} - 1) = {H}")
print("\nFinal Answer String:")
print(f"{m},{n},{H},{P}")