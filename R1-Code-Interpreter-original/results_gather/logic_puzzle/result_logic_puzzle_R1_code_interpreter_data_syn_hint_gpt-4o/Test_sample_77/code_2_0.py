from itertools import combinations

# Given grid
grid = [
    [10, 8, 3, -9, 5],
    [-2, 5, 7, 10, 2],
    [9, 6, 4, -4, 9],
    [4, 4, -6, 4, 1],
    [8, -9, -9, 10, 4]
]

# Start with the given number at (0,4)
initial_position = (0, 4)

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 10 positions from the grid
all_positions = [(r, c) for r in range(5) for c in range(5)]
all_positions.remove(initial_position)

# Try to find a valid combination
for combination in combinations(all_positions, 9):
    selected_positions = [initial_position] + list(combination)
    if is_product_positive(selected_positions, grid):
        print(selected_positions)
        break