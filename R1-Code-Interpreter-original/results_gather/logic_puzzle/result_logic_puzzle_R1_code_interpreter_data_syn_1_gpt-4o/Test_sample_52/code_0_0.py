import itertools

# Define the grid
grid = [
    [1, -4, -1, 5, 2],
    [-7, -6, 7, 4, 2],
    [-2, -5, 4, -1, 6],
    [-3, -4, -7, 4, 4],
    [1, -7, 6, 5, 6]
]

# Start with the given number at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid combination
for combination in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    positions = selected_positions + list(combination)
    if is_valid_selection(positions):
        print(positions)
        break