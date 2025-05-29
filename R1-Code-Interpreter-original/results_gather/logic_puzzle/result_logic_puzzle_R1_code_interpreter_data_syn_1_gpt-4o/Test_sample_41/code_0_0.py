from itertools import combinations

# Define the grid
grid = [
    [-3, -4, 3, 8],
    [2, 3, 1, -2],
    [3, -1, -5, 5],
    [-2, 3, 6, -3]
]

# Pre-select the number 8 at position (0,3)
selected_positions = [(0, 3)]

# Function to calculate the product of selected numbers in a row or column
def product_of_selected(selected, indices, is_row=True):
    product = 1
    for idx in indices:
        if is_row:
            product *= grid[selected[idx][0]][selected[idx][1]]
        else:
            product *= grid[selected[idx][1]][selected[idx][0]]
    return product

# Check if the selection is valid
def is_valid_selection(selected):
    # Check rows
    for row in range(4):
        row_indices = [i for i, (r, c) in enumerate(selected) if r == row]
        if len(row_indices) > 0 and product_of_selected(selected, row_indices) <= 0:
            return False
    # Check columns
    for col in range(4):
        col_indices = [i for i, (r, c) in enumerate(selected) if c == col]
        if len(col_indices) > 0 and product_of_selected(selected, col_indices, is_row=False) <= 0:
            return False
    return True

# Generate all possible combinations of 5 additional positions
all_positions = [(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)]
for additional_positions in combinations(all_positions, 5):
    selected = selected_positions + list(additional_positions)
    if is_valid_selection(selected):
        print(selected)
        break