# Define the grid
grid = [
    [-2, 3, 7, 8, 8],
    [-5, 1, 5, 8, 6],
    [7, -3, -1, -4, 2],
    [5, -4, -7, 3, 3],
    [2, -4, 3, 1, 6]
]

# Selected positions based on analysis
selected_positions = [(0, 4), (0, 1), (0, 2), (1, 1), (1, 2), (2, 0), (2, 4), (3, 3), (4, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Verify the selection
is_valid = is_valid_selection(selected_positions)
print(selected_positions if is_valid else None)