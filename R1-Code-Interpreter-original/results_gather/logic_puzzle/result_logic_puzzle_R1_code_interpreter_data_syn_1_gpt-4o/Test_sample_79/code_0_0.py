import numpy as np
from itertools import combinations

# Define the grid
grid = np.array([
    [-1, 2, -3, 8, -9],
    [2, -8, -5, 1, -3],
    [-10, 6, 1, 1, 10],
    [-2, -2, 6, 9, -6],
    [-3, -10, 2, 6, -9]
])

# Initial selection
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    
    for r, c in selected_positions:
        row_products[r] *= grid[r, c]
        col_products[c] *= grid[r, c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection of 10 numbers
for positions in combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 9):
    selected_positions = [(0, 4)] + list(positions)
    if is_valid_selection(selected_positions):
        break

# Output the selected positions
print(selected_positions)