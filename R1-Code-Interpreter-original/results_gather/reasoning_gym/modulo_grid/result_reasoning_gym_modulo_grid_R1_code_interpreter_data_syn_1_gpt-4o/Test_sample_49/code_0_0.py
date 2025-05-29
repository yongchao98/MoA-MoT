# Positions of the ✅ marks
positions = [
    (1, 15), (2, 16), (3, 5), (4, 9), (5, 4),
    (6, 11), (7, 8), (8, 13), (9, 3), (10, 14),
    (11, 10), (12, 6), (13, 15), (14, 9), (15, 12),
    (16, 2), (17, 3), (18, 17), (19, 15), (20, 17)
]

# Calculate the differences between row and column indices
differences = [row - col for row, col in positions]

# Print the differences
print(differences)