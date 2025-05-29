# Redefine the given checkmark positions
checkmarks = [(1, 3), (6, 16), (8, 4), (12, 13), (13, 14), (18, 11), (20, 3)]

# Calculate modulo of row and column indices with respect to 5
modulo_results = [(row % 5, col % 5) for row, col in checkmarks]

# Print the modulo results to identify a pattern
print(modulo_results)