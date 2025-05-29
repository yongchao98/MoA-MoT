# Define the grid
grid = [
    [-4, -5, -7, -6, 6],
    [5, 7, 4, -6, -1],
    [2, -5, -7, 5, 6],
    [4, -6, -3, 3, -3],
    [9, 5, -5, 5, 6]
]

# Manually selected positions
selected_positions = [(0, 4), (1, 0), (1, 1), (2, 0), (2, 3), (3, 0), (3, 3), (4, 0), (4, 1), (0, 0)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Check if the selection meets the condition
if is_product_positive(selected_positions, grid):
    # Format the output
    result = ', '.join(f"{r} {c}" for r, c in selected_positions)
    print(f"<<<{result}>>>")
else:
    print("No valid selection found")