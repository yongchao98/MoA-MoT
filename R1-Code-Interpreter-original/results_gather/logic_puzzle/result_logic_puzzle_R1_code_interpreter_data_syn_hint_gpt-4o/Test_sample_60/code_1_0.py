# Given grid
grid = [
    [7, 4, -2, -8, 8],
    [2, 2, 4, 4, 1],
    [4, -5, 5, -1, 3],
    [-5, 3, 6, 4, 2],
    [6, 6, -4, -6, 2]
]

# Selected positions from previous code
selected_positions = [(0, 4), (0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (4, 2), (4, 3)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Verify the selection
if is_valid_selection(selected_positions):
    # Format the output
    formatted_output = ', '.join(f"{r} {c}" for r, c in selected_positions)
    print(f"<<<{formatted_output}>>>")
else:
    print("No valid selection found.")