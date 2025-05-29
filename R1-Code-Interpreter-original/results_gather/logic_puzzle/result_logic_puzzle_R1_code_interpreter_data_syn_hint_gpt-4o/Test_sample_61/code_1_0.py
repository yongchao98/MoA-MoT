# Define the grid
grid = [
    [-8, 8, -4, 1, -6],
    [-7, 7, -1, 5, 4],
    [3, -1, -3, 3, -2],
    [2, 4, 8, 8, -3],
    [-6, -7, -7, -6, 3]
]

# Start with the given number at (0,4)
selected_positions = [(0, 4)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    # Check row and column products
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection
def find_selection():
    # Try all combinations of 9 positions
    from itertools import combinations
    all_positions = [(r, c) for r in range(5) for c in range(5)]
    for positions in combinations(all_positions, 9):
        if (0, 4) in positions:
            if is_valid_selection(positions):
                return positions
    return None

# Find a valid selection
valid_selection = find_selection()

# Format the output
if valid_selection:
    result = ', '.join(f"{r} {c}" for r, c in valid_selection)
else:
    result = "No valid selection found"

print(result)