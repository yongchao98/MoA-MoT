import itertools

grid = [
    [-13, 6, 12, 14, 12, -3, -7],
    [13, -13, 1, 13, 13, 2, -10],
    [2, 6, -9, 13, 1, 10, 12],
    [8, -6, 13, 8, 14, 9, 13],
    [3, 10, -14, 10, -6, 5, 11],
    [-9, 11, 4, 5, 7, -11, -10],
    [-5, 1, 4, 13, 14, 2, 12]
]

# Start with the given number
initial_position = (0, 6)
selected_positions = [initial_position]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_products = [1] * 7
    col_products = [1] * 7
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Backtracking function to find a valid selection
def find_selection(selected_positions, grid, start_row=0, start_col=0):
    if len(selected_positions) == 16:
        if is_product_positive(selected_positions, grid):
            return selected_positions
        return None

    for r in range(start_row, 7):
        for c in range(start_col if r == start_row else 0, 7):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                result = find_selection(selected_positions, grid, r, c + 1)
                if result:
                    return result
                selected_positions.pop()
    return None

# Find the selection
result = find_selection(selected_positions, grid)

# Print the selected positions
print(result)