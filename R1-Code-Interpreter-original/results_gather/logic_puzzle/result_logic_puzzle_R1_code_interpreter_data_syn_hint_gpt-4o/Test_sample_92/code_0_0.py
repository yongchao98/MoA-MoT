# Define the grid
grid = [
    [-4, 2, 3, -10, 11, 2],
    [5, 4, -6, -10, 4, 4],
    [6, -9, 4, 12, 6, 8],
    [12, -10, 6, 6, -5, 1],
    [11, -4, 1, 11, 11, 7],
    [-5, -3, 10, 11, -6, 12]
]

# Start with the given number at position (0,5)
selected_positions = [(0, 5)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    # Check row and column products
    row_products = [1] * 6
    col_products = [1] * 6
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to select 12 numbers
def select_numbers():
    # Start with the initial selection
    selected_positions = [(0, 5)]
    # Iterate over the grid to select numbers
    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if len(selected_positions) == 12:
                    if is_valid_selection(selected_positions):
                        return selected_positions
                    selected_positions.pop()
    return None

# Get the selected positions
selected_positions = select_numbers()

# Print the selected positions
print(selected_positions)