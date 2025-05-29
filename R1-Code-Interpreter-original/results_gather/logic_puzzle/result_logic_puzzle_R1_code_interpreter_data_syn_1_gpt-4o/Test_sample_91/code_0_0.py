import itertools

# Define the grid
grid = [
    [12, 3, 5, 8, 1, 8],
    [12, 8, 1, 3, 1, 2],
    [5, 4, -9, -8, 9, -4],
    [7, 10, -6, -7, 7, 6],
    [3, 1, 10, 9, -10, 12],
    [2, 11, -4, 8, 11, -12]
]

# Pre-select the number at (0,5)
selected_positions = [(0, 5)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 6
    col_products = [1] * 6
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection of 12 numbers
for positions in itertools.combinations([(r, c) for r in range(6) for c in range(6) if (r, c) != (0, 5)], 11):
    current_selection = selected_positions + list(positions)
    if is_valid_selection(current_selection):
        print(current_selection)
        break