# Define the grid
grid = [
    [7, 8, 1, 9, 4],
    [-6, 1, 6, 5, 6],
    [-7, 8, 5, -9, 6],
    [10, 1, 2, 3, 9],
    [5, 3, -10, 6, 10]
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

# Function to count negative numbers in a row or column
def count_negatives(selected_positions, grid):
    row_negatives = [0] * 5
    col_negatives = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_negatives[r] += 1
            col_negatives[c] += 1
    return row_negatives, col_negatives

# Try to select 9 more numbers
for r in range(5):
    for c in range(5):
        if (r, c) != (0, 4) and (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 10:
                row_negatives, col_negatives = count_negatives(selected_positions, grid)
                if all(n % 2 == 0 for n in row_negatives) and all(n % 2 == 0 for n in col_negatives):
                    if is_product_positive(selected_positions, grid):
                        break
                selected_positions.pop()
    if len(selected_positions) == 10:
        break

# Output the selected positions
print(selected_positions)