import itertools

# Define the grid
grid = [
    [7, -3, 5, 6],
    [6, 6, -3, 2],
    [4, 4, 5, -4],
    [2, 8, 8, 7]
]

# Pre-select the number at (0,3)
selected_positions = [(0, 3)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(positions):
    row_products = [1] * 4
    col_products = [1] * 4
    
    for r, c in positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products if p != 1) and all(p > 0 for p in col_products if p != 1)

# Try to find a valid selection
for positions in itertools.combinations([(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)], 5):
    positions = selected_positions + list(positions)
    if is_valid_selection(positions):
        print(positions)
        break