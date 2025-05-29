# Given positions of "✅" marks
positions = [
    (1, 8), (2, 13), (3, 16), (4, 17), (5, 18),
    (6, 10), (7, 19), (8, 2), (9, 9), (10, 5),
    (11, 14), (12, 11), (13, 19), (14, 3), (15, 10),
    (16, 4), (17, 5), (18, 7), (19, 13)
]

# Extract the column positions
columns = [pos[1] for pos in positions]

# Calculate the differences between consecutive column positions
differences = [columns[i+1] - columns[i] for i in range(len(columns) - 1)]

# Print the differences to analyze the pattern
print("Differences:", differences)