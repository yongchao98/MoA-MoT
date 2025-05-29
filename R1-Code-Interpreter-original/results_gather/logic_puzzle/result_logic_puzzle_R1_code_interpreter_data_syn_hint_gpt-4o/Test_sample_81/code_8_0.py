grid = [
    [-2, 9, -9, -3, -3, 11],
    [10, 10, 9, 6, 5, -11],
    [5, 11, 9, 9, 1, 8],
    [4, 3, -5, -11, 4, 10],
    [11, -8, 6, 6, 2, 7],
    [-12, -11, 6, 2, 5, -3]
]

# Start with the given number
selected_positions = [(0, 5)]
selected_numbers = [11]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_products = [1] * 6
    col_products = [1] * 6
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Backtracking function to find a valid selection
def find_selection(selected_positions, grid):
    if len(selected_positions) == 12:
        if is_product_positive(selected_positions, grid):
            return selected_positions
        else:
            return None

    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if is_product_positive(selected_positions, grid):
                    result = find_selection(selected_positions, grid)
                    if result is not None:
                        return result
                selected_positions.pop()
    return None

# Find the selection
result = find_selection(selected_positions, grid)

# Output the selected positions
print(result)