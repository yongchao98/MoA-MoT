# Define the grid
grid = [
    [12, 3, 5, 8, 1, 8],
    [12, 8, 1, 3, 1, 2],
    [5, 4, -9, -8, 9, -4],
    [7, 10, -6, -7, 7, 6],
    [3, 1, 10, 9, -10, 12],
    [2, 11, -4, 8, 11, -12]
]

# Start with the given number
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_count = [0] * 6
    col_count = [0] * 6
    row_product = [1] * 6
    col_product = [1] * 6
    
    for r, c in selected_positions:
        row_count[r] += 1
        col_count[c] += 1
        row_product[r] *= grid[r][c]
        col_product[c] *= grid[r][c]
    
    for i in range(6):
        if row_count[i] > 0 and row_product[i] <= 0:
            return False
        if col_count[i] > 0 and col_product[i] <= 0:
            return False
    return True

# Select numbers ensuring the product condition
for r in range(6):
    for c in range(6):
        if len(selected_positions) >= 12:
            break
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if not is_product_positive(selected_positions, grid):
                selected_positions.pop()

# Print the selected positions
print(selected_positions)