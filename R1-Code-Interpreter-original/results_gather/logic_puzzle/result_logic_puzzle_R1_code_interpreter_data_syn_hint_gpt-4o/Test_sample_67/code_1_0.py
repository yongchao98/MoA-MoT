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

# Backtracking function to find valid selections
def find_selections(selected_positions, grid, start_row, start_col):
    if len(selected_positions) == 10:
        if is_product_positive(selected_positions, grid):
            return selected_positions
        else:
            return None

    for r in range(start_row, 5):
        for c in range(start_col if r == start_row else 0, 5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                result = find_selections(selected_positions, grid, r, c + 1)
                if result:
                    return result
                selected_positions.pop()
    return None

# Find the valid selection
result = find_selections(selected_positions, grid, 0, 0)

# Output the selected positions
print(result)