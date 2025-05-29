# Define the positions of the "✅" marks
positions = [
    (1, 15), (2, 16), (3, 5), (4, 9), (5, 4),
    (6, 12), (7, 8), (8, 7), (9, 4), (10, 13),
    (11, 10), (12, 6), (13, 14), (14, 9), (15, 12),
    (16, 2), (17, 3), (18, 18), (19, 15), (20, 17)
]

# Extract column indices
cols = [col for row, col in positions]

# Calculate the differences between consecutive column indices
column_differences = [cols[i] - cols[i - 1] for i in range(1, len(cols))]

# Print the differences
print("Column Differences:", column_differences)