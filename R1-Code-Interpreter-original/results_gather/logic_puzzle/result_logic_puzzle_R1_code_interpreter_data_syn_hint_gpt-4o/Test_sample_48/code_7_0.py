import itertools

# Define the grid
grid = [
    [-7, 1, 1, 6, 5],
    [7, 1, 3, 8, 8],
    [7, -1, -2, -5, 1],
    [7, -1, 6, -7, 2],
    [3, -3, 5, -8, 7]
]

# Pre-select the number 5 at position (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 8 additional positions
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)]
for additional_positions in itertools.combinations(all_positions, 8):
    current_selection = selected_positions + list(additional_positions)
    if is_valid_selection(current_selection):
        print(current_selection)
        break