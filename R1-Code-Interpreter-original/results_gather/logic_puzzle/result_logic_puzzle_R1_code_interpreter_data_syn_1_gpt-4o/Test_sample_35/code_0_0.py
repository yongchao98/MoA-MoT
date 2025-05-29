# Define the grid
grid = [
    [1, 3, -3, 8],
    [8, 7, 8, -3],
    [2, 2, -4, -1],
    [8, 2, 2, 4]
]

# Start with the given number at position (0,3)
selected_positions = [(0, 3)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    # Check row and column products
    row_products = [1] * 4
    col_products = [1] * 4
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to select 5 more numbers
for r in range(4):
    for c in range(4):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 6 and is_valid_selection(selected_positions):
                print(selected_positions)
                break
            if len(selected_positions) < 6:
                continue
            selected_positions.pop()
    if len(selected_positions) == 6:
        break