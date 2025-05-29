# Define the positions of the checkmarks
checkmark_positions = [
    (1, 7), (2, 12), (3, 8), (4, 6), (5, 15), (6, 4), (7, 1), (7, 18),
    (8, 3), (9, 14), (10, 16), (11, 13), (12, 2), (12, 19), (13, 11),
    (14, 9), (15, 5), (16, 10), (19, 7), (20, 13)
]

# Check for a pattern in the positions
def find_pattern(positions):
    # Check for a simple arithmetic pattern
    for i in range(1, len(positions)):
        row_diff = positions[i][0] - positions[i-1][0]
        col_diff = positions[i][1] - positions[i-1][1]
        print(f"Row difference: {row_diff}, Column difference: {col_diff}")

find_pattern(checkmark_positions)