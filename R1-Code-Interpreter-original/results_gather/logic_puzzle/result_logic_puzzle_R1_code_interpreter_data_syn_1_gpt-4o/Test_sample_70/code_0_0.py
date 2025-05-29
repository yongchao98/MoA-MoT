import itertools

# Define the grid
grid = [
    [-2, -3, -8, 8, 9],
    [-5, 7, -3, 3, 2],
    [-6, -2, 4, 6, -8],
    [8, -10, 10, 8, -6],
    [5, -9, 8, -10, -4]
]

# Initial selected position
initial_position = (0, 4)

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 9 additional positions
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != initial_position]
combinations = itertools.combinations(all_positions, 9)

# Find a valid combination
for combo in combinations:
    selected_positions = [initial_position] + list(combo)
    if is_valid_selection(selected_positions):
        result = selected_positions
        break

# Format the result
formatted_result = ', '.join(f"{r} {c}" for r, c in result)
print(formatted_result)