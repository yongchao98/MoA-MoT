# Define the positions of the "✅" marks
positions = [(1, 17), (6, 8), (8, 6), (12, 13), (13, 11), (14, 12), (18, 2), (18, 19), (19, 17)]

# Calculate the differences between row and column indices
differences = [(row - col) for row, col in positions]

# Print the differences
print(differences)