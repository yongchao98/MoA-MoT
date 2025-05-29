import itertools

# Define the grid
grid = [
    [2, -10, -3, 1, -1],
    [5, 6, 6, -9, 3],
    [9, -4, -7, 10, -5],
    [10, -2, 4, 1, -10],
    [-8, 6, -2, 10, 1]
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
    if is_valid_selection(selected_positions + list(positions)):
        selected_positions += list(positions)
        break

# Print the selected positions
print(selected_positions)