# Calculate modulo of row and column indices
modulo_results = [(row % 5, col % 5) for row, col in checkmarks]

# Print the modulo results to identify a pattern
print(modulo_results)