# Positions of '✅' in the grid
positions = [
    (1, 18), (2, 9), (2, 19), (3, 6), (6, 3), (6, 13),
    (7, 14), (9, 2), (11, 18), (13, 6), (14, 7), (14, 17),
    (18, 14), (19, 1), (19, 11)
]

# Analyze the differences between row and column indices
differences = [(r - c) for r, c in positions]

# Check if the position (1, 4) fits into the pattern
def check_position(row, col, differences):
    difference = row - col
    return difference in differences

# Check the position of the '❔'
position_check = check_position(1, 4, differences)

# Print the result of the check
print(position_check)