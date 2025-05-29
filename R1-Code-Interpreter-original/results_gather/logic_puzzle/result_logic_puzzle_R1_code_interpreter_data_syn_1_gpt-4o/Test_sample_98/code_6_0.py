# Define the grid
grid = [
    [-13, 6, 12, 14, 12, -3, -7],
    [13, -13, 1, 13, 13, 2, -10],
    [2, 6, -9, 13, 1, 10, 12],
    [8, -6, 13, 8, 14, 9, 13],
    [3, 10, -14, 10, -6, 5, 11],
    [-9, 11, 4, 5, 7, -11, -10],
    [-5, 1, 4, 13, 14, 2, 12]
]

# Start with the given number -7 at position (0,6)
selected_positions = [(0, 6)]

# Manually select numbers to ensure even negative count in each row and column
# This is a manual selection based on trial and error to meet the criteria
manual_selection = [
    (0, 6), (0, 1), (0, 3), (0, 4),
    (1, 0), (1, 3), (1, 4), (1, 5),
    (2, 0), (2, 1), (2, 3), (2, 6),
    (3, 0), (3, 2), (3, 3), (3, 6)
]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 7
    col_neg_count = [0] * 7
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    # Check if all rows and columns have even negative counts
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Check if the manual selection is valid
is_valid = is_valid_selection(manual_selection)

# Print the selected positions if valid
if is_valid:
    print(manual_selection)
else:
    print("No valid selection found")