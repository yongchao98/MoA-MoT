# Define the grid
grid = [
    [-9, 11, 7, 9, 8, 5],
    [1, 3, -3, -3, -4, 12],
    [1, 7, 6, 11, 10, 2],
    [7, 9, 3, 1, 7, 8],
    [11, -3, 5, 12, 1, 9],
    [-3, 8, -1, 6, -5, 5]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 6
    col_products = [1] * 6
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Greedily select numbers to maintain positive products
for r in range(6):
    for c in range(6):
        if (r, c) != (0, 5) and len(selected_positions) < 12:
            selected_positions.append((r, c))
            if not is_valid_selection(selected_positions):
                selected_positions.pop()

# Print the selected positions
print(selected_positions)