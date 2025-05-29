import itertools

# Define the grid
grid = [
    [-4, 8, 4, 4],
    [4, -5, 4, -2],
    [-1, 6, 2, 8],
    [5, -3, 6, 3]
]

# Initial selected position
initial_position = (0, 3)

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    # Initialize row and column products
    row_products = [1] * 4
    col_products = [1] * 4
    
    # Calculate the product for each row and column
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    # Check if all products are positive
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 5 additional positions
all_positions = [(r, c) for r in range(4) for c in range(4) if (r, c) != initial_position]
combinations = itertools.combinations(all_positions, 5)

# Find a valid combination
for combo in combinations:
    selected_positions = [initial_position] + list(combo)
    if is_valid_selection(selected_positions):
        # Format the output as required
        result = ', '.join(f"{r} {c}" for r, c in selected_positions)
        print(result)
        break