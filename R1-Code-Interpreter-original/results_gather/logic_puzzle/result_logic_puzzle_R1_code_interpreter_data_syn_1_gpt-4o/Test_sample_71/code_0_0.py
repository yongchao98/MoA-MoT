import itertools

# Define the grid
grid = [
    [-8, 5, -9, 6, 7],
    [10, -7, 3, 9, 9],
    [-2, 7, -5, 10, 1],
    [4, 9, 6, -5, 1],
    [7, 9, 3, 7, 5]
]

# Initial selected position
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection of 10 numbers
def find_selection():
    # Generate all possible positions in the grid
    all_positions = [(r, c) for r in range(5) for c in range(5)]
    # Remove the initial selected position
    all_positions.remove((0, 4))
    
    # Try all combinations of 9 additional positions
    for additional_positions in itertools.combinations(all_positions, 9):
        current_selection = selected_positions + list(additional_positions)
        if is_valid_selection(current_selection):
            return current_selection
    return None

# Find a valid selection
valid_selection = find_selection()

# Format the output
if valid_selection:
    result = ', '.join(f"{r} {c}" for r, c in valid_selection)
    print(f"<<<{result}>>>")
else:
    print("No valid selection found.")