# Number of holes, which corresponds to the genus (g) of the surface.
g = 2

# Number of separate patches, which corresponds to the number of boundaries (b).
b = 2

# The formula for the number of logical qubits (k) is k = 2g + b - 1.
k = 2 * g + b - 1

# Print the final equation with the values plugged in.
print(f"The number of logical qubits is calculated as: 2 * {g} + {b} - 1 = {k}")