import itertools

# Given grid
grid = [
    [7, 4, -2, -8, 8],
    [2, 2, 4, 4, 1],
    [4, -5, 5, -1, 3],
    [-5, 3, 6, 4, 2],
    [6, 6, -4, -6, 2]
]

# Start with the given number 8 at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection of 9 numbers
for positions in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    positions = selected_positions + list(positions)
    if is_valid_selection(positions):
        print(positions)
        break