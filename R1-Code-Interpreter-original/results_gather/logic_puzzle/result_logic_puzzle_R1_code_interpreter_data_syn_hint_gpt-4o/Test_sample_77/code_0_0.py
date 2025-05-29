# Given grid
grid = [
    [10, 8, 3, -9, 5],
    [-2, 5, 7, 10, 2],
    [9, 6, 4, -4, 9],
    [4, 4, -6, 4, 1],
    [8, -9, -9, 10, 4]
]

# Start with the given number at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to select 10 numbers
for r in range(5):
    for c in range(5):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 10:
                if is_product_positive(selected_positions, grid):
                    break
                else:
                    selected_positions.pop()
    if len(selected_positions) == 10:
        break

# Output the selected positions
print(selected_positions)