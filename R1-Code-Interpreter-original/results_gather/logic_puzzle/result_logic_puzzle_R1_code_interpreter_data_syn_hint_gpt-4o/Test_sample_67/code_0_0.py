# Given grid
grid = [
    [4, 7, -7, 6, 2],
    [-6, -3, 3, 10, 2],
    [9, 4, 10, 5, 4],
    [-4, 4, 8, 5, 7],
    [8, 1, 10, -2, 10]
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
    return all(p > 0 for p in row_products if p != 1) and all(p > 0 for p in col_products if p != 1)

# Select numbers ensuring the product condition
for r in range(5):
    for c in range(5):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 10:
                if is_product_positive(selected_positions, grid):
                    break
                else:
                    selected_positions.pop()
            elif not is_product_positive(selected_positions, grid):
                selected_positions.pop()
    if len(selected_positions) == 10:
        break

# Output the selected positions
print(selected_positions)