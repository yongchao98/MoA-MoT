import itertools

# Define the grid
grid = [
    [12, 2, 8, 9, 4, 10],
    [-9, -12, 10, 2, -8, 7],
    [-4, 9, 4, -8, -10, -11],
    [-4, 12, -1, 11, -5, -2],
    [-1, 8, 11, 3, -6, 2],
    [11, -6, -12, -11, 4, 4]
]

# Initial selected position
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
def find_selection():
    # Generate all possible positions except the initial one
    all_positions = [(r, c) for r in range(6) for c in range(6) if (r, c) != (0, 5)]
    
    # Try all combinations of 11 positions from the remaining positions
    for combination in itertools.combinations(all_positions, 11):
        current_selection = selected_positions + list(combination)
        if is_valid_selection(current_selection):
            return current_selection
    return None

# Find the selection
selection = find_selection()

# Format the output
if selection:
    result = ', '.join(f"{r} {c}" for r, c in selection)
    print(f"<<<{result}>>>")
else:
    print("No valid selection found.")