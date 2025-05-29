# Define the positions of the "✅" symbols
positions = [
    (1, 14), (2, 8), (4, 11), (8, 3), (8, 17), 
    (9, 13), (12, 5), (12, 19), (13, 9), (14, 8), 
    (15, 2), (15, 16), (17, 14), (18, 8), (20, 12)
]

# Check for a pattern in the differences between consecutive positions
def find_differences(positions):
    row_diffs = []
    col_diffs = []
    for i in range(1, len(positions)):
        row_diff = positions[i][0] - positions[i-1][0]
        col_diff = positions[i][1] - positions[i-1][1]
        row_diffs.append(row_diff)
        col_diffs.append(col_diff)
    return row_diffs, col_diffs

row_diffs, col_diffs = find_differences(positions)
print("Row differences:", row_diffs)
print("Column differences:", col_diffs)