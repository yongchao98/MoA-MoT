# Define the grid
grid = [
    [2, -2, -5, -7, -2],
    [7, 2, -7, 8, -1],
    [4, -5, -6, 5, 6],
    [7, -6, 2, 8, 4],
    [-1, 8, 4, -4, 5]
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

# Try to select 9 numbers
def select_numbers():
    # Start with the given position
    selected_positions = [(0, 4)]
    # Try to add more positions
    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if len(selected_positions) == 9 and is_valid_selection(selected_positions):
                    return selected_positions
                if len(selected_positions) < 9:
                    continue
                selected_positions.pop()
    return None

# Get the selected positions
selected_positions = select_numbers()

# Print the result
print(selected_positions)