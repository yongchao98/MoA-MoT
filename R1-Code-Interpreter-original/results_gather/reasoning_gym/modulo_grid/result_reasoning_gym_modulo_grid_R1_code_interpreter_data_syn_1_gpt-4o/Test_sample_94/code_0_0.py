# Given checkmark positions
checkmarks = [(1, 3), (6, 16), (8, 4), (12, 13), (13, 14), (18, 11), (20, 3)]

# Calculate differences between row and column indices
differences = [row - col for row, col in checkmarks]

# Print the differences to identify a pattern
print(differences)