from itertools import combinations

# Define the grid
grid = [
    [8, 3, -8, -6, 5],
    [-2, -7, 7, -3, 3],
    [3, 2, 2, 8, -3],
    [4, 5, -5, -2, -5],
    [4, -4, -7, 5, 8]
]

# Pre-select the given number
pre_selected = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 8 additional positions
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)]
for additional_positions in combinations(all_positions, 8):
    selected_positions = pre_selected + list(additional_positions)
    if is_valid_selection(selected_positions):
        # Format the output as required
        result = ', '.join(f"{r} {c}" for r, c in selected_positions)
        print(f"<<<{result}>>>")
        break