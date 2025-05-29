from itertools import combinations

# Define the grid
grid = [
    [-5, 1, -3, 7, 5],
    [6, 5, 2, 8, 3],
    [7, 8, 8, -7, 6],
    [-4, 2, -6, -7, -4],
    [-3, 5, -4, 8, 3]
]

# Pre-select the given number
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 8 more positions
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)]
for combo in combinations(all_positions, 8):
    current_selection = selected_positions + list(combo)
    if is_valid_selection(current_selection):
        print(current_selection)
        break