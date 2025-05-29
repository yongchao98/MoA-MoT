from itertools import combinations

# Define the grid
grid = [
    [-2, 3, 7, 8, 8],
    [-5, 1, 5, 8, 6],
    [7, -3, -1, -4, 2],
    [5, -4, -7, 3, 3],
    [2, -4, 3, 1, 6]
]

# Predefined selected position
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection
for positions in combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    current_selection = selected_positions + list(positions)
    if is_valid_selection(current_selection):
        print(current_selection)
        break