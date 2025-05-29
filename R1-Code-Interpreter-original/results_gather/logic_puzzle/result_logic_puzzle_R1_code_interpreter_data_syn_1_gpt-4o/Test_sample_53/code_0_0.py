import itertools

# The given grid
grid = [
    [3, 4, -7, 7, 3],
    [-5, -2, 2, -1, 1],
    [-6, 4, 1, -4, 7],
    [7, 6, 1, 3, -6],
    [6, -3, 3, -3, 8]
]

# The initial selected position
initial_position = (0, 4)

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 9 positions from the grid
all_positions = [(r, c) for r in range(5) for c in range(5)]
all_positions.remove(initial_position)  # Remove the initial position from the list

# Try all combinations of 8 positions (since we already have 1 position)
for combination in itertools.combinations(all_positions, 8):
    selected_positions = [initial_position] + list(combination)
    if is_valid_selection(selected_positions):
        print(selected_positions)
        break