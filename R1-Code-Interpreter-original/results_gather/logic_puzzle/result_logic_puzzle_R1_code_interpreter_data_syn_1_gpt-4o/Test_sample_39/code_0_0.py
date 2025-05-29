import itertools

# Define the grid
grid = [
    [8, 6, 2, -1],
    [1, 8, -3, 5],
    [4, -4, 1, 4],
    [8, 7, -5, 3]
]

# Pre-select the number at (0,3)
selected_positions = [(0, 3)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 4
    col_products = [1] * 4
    
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try all combinations of 5 more positions
all_positions = [(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)]
for combination in itertools.combinations(all_positions, 5):
    current_selection = selected_positions + list(combination)
    if is_valid_selection(current_selection):
        print(current_selection)
        break