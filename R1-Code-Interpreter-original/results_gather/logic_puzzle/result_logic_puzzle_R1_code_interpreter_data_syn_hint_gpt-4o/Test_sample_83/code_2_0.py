grid = [
    [8, -2, 7, 6, -4, 3],
    [11, -11, 7, 11, -8, 8],
    [-4, -6, 10, -6, 4, -12],
    [3, 10, 1, 9, -1, -4],
    [4, 7, 12, -10, 7, 11],
    [8, 7, 9, -11, 3, 7]
]

# Start with the given number at (0,5)
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_products = [1] * 6
    col_products = [1] * 6
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Function to count negative numbers in a row or column
def count_negatives(selected_positions, grid):
    row_negatives = [0] * 6
    col_negatives = [0] * 6
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_negatives[r] += 1
            col_negatives[c] += 1
    return row_negatives, col_negatives

# Try to select 12 numbers
def select_numbers(grid, selected_positions):
    if len(selected_positions) == 12:
        if is_product_positive(selected_positions, grid):
            return selected_positions
        return None

    row_negatives, col_negatives = count_negatives(selected_positions, grid)

    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                # Check if adding this number keeps the negatives even
                if (row_negatives[r] % 2 == 0 and grid[r][c] < 0) or (row_negatives[r] % 2 == 1 and grid[r][c] > 0):
                    if (col_negatives[c] % 2 == 0 and grid[r][c] < 0) or (col_negatives[c] % 2 == 1 and grid[r][c] > 0):
                        selected_positions.append((r, c))
                        result = select_numbers(grid, selected_positions)
                        if result:
                            return result
                        selected_positions.pop()
    return None

# Find the solution
solution = select_numbers(grid, selected_positions)
print(solution)