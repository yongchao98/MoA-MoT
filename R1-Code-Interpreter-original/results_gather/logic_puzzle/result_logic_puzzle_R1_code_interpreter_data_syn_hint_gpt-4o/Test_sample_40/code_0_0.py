from itertools import combinations

# Define the grid
grid = [
    [-4, 5, 7, -1],
    [2, -4, -2, -5],
    [6, 6, -4, 4],
    [-1, 2, 1, -1]
]

# Pre-selected position
selected_positions = [(0, 3)]

# Function to calculate product of selected numbers in a row or column
def product_of_selected(selected, grid, is_row=True):
    products = [1] * 4
    for r, c in selected:
        if is_row:
            products[r] *= grid[r][c]
        else:
            products[c] *= grid[r][c]
    return products

# Check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected, grid):
    row_products = product_of_selected(selected, grid, is_row=True)
    col_products = product_of_selected(selected, grid, is_row=False)
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 5 additional positions
all_positions = [(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)]
for additional_positions in combinations(all_positions, 5):
    selected = selected_positions + list(additional_positions)
    if is_valid_selection(selected, grid):
        # Format the output as required
        result = ', '.join(f"{r} {c}" for r, c in selected)
        print(f"<<<{result}>>>")
        break