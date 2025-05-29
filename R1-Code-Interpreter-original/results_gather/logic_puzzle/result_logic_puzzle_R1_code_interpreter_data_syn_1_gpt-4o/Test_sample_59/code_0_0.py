from itertools import combinations

# Define the grid
grid = [
    [5, -2, 2, -8, -6],
    [1, -2, 7, -1, 8],
    [4, -7, 6, -7, -3],
    [8, -7, 6, 8, 8],
    [-4, -5, 4, 6, 3]
]

# Predefined selection
preselected = [(0, 4)]

# Function to check if a selection is valid
def is_valid_selection(selection):
    # Check row and column products
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selection:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    # Check if all products are positive
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible selections of 9 numbers including the preselected one
all_positions = [(r, c) for r in range(5) for c in range(5)]
all_positions.remove((0, 4))  # Remove the preselected position

# Try all combinations
for combo in combinations(all_positions, 8):
    selection = preselected + list(combo)
    if is_valid_selection(selection):
        print(selection)
        break