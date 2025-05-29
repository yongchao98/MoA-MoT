# Define the positions of the "✅" marks
positions = [(1, 17), (6, 8), (8, 6), (12, 13), (13, 11), (14, 12), (18, 2), (18, 19), (19, 17)]

# Calculate the modulo of row and column indices with 5
modulo_results = [(row % 5, col % 5) for row, col in positions]

# Print the modulo results
print(modulo_results)