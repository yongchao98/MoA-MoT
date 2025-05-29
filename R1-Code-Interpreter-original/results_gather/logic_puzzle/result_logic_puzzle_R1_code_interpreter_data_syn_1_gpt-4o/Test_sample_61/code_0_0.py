import itertools

# Define the grid
grid = [
    [-8,  8, -4,  1, -6],
    [-7,  7, -1,  5,  4],
    [ 3, -1, -3,  3, -2],
    [ 2,  4,  8,  8, -3],
    [-6, -7, -7, -6,  3]
]

# Pre-select the first number
selected_positions = [(0, 4)]

# Function to calculate the product of selected numbers in a row or column
def product_of_selected(selected, indices, is_row=True):
    prod = 1
    for idx in indices:
        if is_row:
            prod *= grid[selected[idx][0]][selected[idx][1]]
        else:
            prod *= grid[selected[idx][1]][selected[idx][0]]
    return prod

# Check if the selection is valid
def is_valid_selection(selected):
    # Check rows
    for row in range(5):
        row_indices = [i for i, (r, c) in enumerate(selected) if r == row]
        if row_indices and product_of_selected(selected, row_indices, is_row=True) <= 0:
            return False
    # Check columns
    for col in range(5):
        col_indices = [i for i, (r, c) in enumerate(selected) if c == col]
        if col_indices and product_of_selected(selected, col_indices, is_row=False) <= 0:
            return False
    return True

# Try to find a valid selection
for combination in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    selected = selected_positions + list(combination)
    if is_valid_selection(selected):
        print(selected)
        break