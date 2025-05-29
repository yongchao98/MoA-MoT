from itertools import combinations

# Define the grid
grid = [
    [2, -2, -5, -7, -2],
    [7, 2, -7, 8, -1],
    [4, -5, -6, 5, 6],
    [7, -6, 2, 8, 4],
    [-1, 8, 4, -4, 5]
]

# Start with the given number at (0,4)
initial_position = (0, 4)

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    # Check row and column products
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 9 positions including the initial position
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != initial_position]
for combination in combinations(all_positions, 8):
    selected_positions = [initial_position] + list(combination)
    if is_valid_selection(selected_positions):
        print(selected_positions)
        break