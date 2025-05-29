from itertools import combinations

# Define the grid
grid = [
    [1, 2, 4, -2],
    [-4, 3, 3, 1],
    [2, 1, 7, 2],
    [7, -3, 7, -3]
]

# Pre-select the given number
selected_positions = [(0, 3)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 4
    col_products = [1] * 4
    
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Generate all possible combinations of 5 additional positions
all_positions = [(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)]
for additional_positions in combinations(all_positions, 5):
    current_selection = selected_positions + list(additional_positions)
    if is_valid_selection(current_selection):
        # Convert positions to the required format
        result = ', '.join(f"{r} {c}" for r, c in current_selection)
        print(f"<<<{result}>>>")
        break