# Redefine the given checkmark positions
checkmarks = [(1, 3), (6, 16), (8, 4), (12, 13), (13, 14), (18, 11), (20, 3)]

# Calculate sums of row and column indices
sums = [row + col for row, col in checkmarks]

# Print the sums to identify a pattern
print(sums)