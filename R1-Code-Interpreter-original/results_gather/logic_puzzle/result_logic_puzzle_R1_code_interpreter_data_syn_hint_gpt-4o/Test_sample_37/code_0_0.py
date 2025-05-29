from itertools import combinations

# Define the grid
grid = [
    [-1, 8, 4, 1],
    [3, 4, 6, 1],
    [-4, 1, -1, 4],
    [-1, 2, 4, -2]
]

# Pre-select the number at (0, 3)
selected_positions = [(0, 3)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 4
    col_products = [1] * 4
    
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection
for positions in combinations([(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)], 5):
    current_selection = selected_positions + list(positions)
    if is_valid_selection(current_selection):
        print(current_selection)
        break