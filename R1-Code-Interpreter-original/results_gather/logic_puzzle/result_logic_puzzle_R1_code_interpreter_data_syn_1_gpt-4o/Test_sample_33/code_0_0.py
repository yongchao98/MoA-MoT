import itertools

# Define the grid
grid = [
    [6, 1, -5, -5],
    [8, -4, -4, 2],
    [-4, 2, 8, 5],
    [4, 8, -2, 5]
]

# Initial selected position
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
for positions in itertools.combinations([(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)], 5):
    selected_positions = [(0, 3)] + list(positions)
    if is_valid_selection(selected_positions):
        break

# Print the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))