import itertools

# Define the grid
grid = [
    [-9, 11, 7, 9, 8, 5],
    [1, 3, -3, -3, -4, 12],
    [1, 7, 6, 11, 10, 2],
    [7, 9, 3, 1, 7, 8],
    [11, -3, 5, 12, 1, 9],
    [-3, 8, -1, 6, -5, 5]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 6
    col_products = [1] * 6
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to select 12 numbers
for positions in itertools.combinations([(r, c) for r in range(6) for c in range(6) if (r, c) != (0, 5)], 11):
    current_selection = selected_positions + list(positions)
    if is_valid_selection(current_selection):
        selected_positions = current_selection
        break

# Print the selected positions
print(selected_positions)