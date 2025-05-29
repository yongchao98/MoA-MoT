import itertools

# Define the grid
grid = [
    [7, 8, 1, 9, 4],
    [-6, 1, 6, 5, 6],
    [-7, 8, 5, -9, 6],
    [10, 1, 2, 3, 9],
    [5, 3, -10, 6, 10]
]

# Pre-select the number at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection
for positions in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 9):
    current_selection = selected_positions + list(positions)
    if is_valid_selection(current_selection):
        print(current_selection)
        break