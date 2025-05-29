# Define the grid
grid = [
    [-2, 9, -9, -3, -3, 11],
    [10, 10, 9, 6, 5, -11],
    [5, 11, 9, 9, 1, 8],
    [4, 3, -5, -11, 4, 10],
    [11, -8, 6, 6, 2, 7],
    [-12, -11, 6, 2, 5, -3]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_count = [0] * 6
    col_count = [0] * 6
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    # Check if all rows and columns have an even number of negative numbers
    return all(count % 2 == 0 for count in row_count) and all(count % 2 == 0 for count in col_count)

# Try to select 12 numbers
for r in range(6):
    for c in range(6):
        if (r, c) != (0, 5) and len(selected_positions) < 12:
            selected_positions.append((r, c))
            if not is_product_positive(selected_positions, grid):
                selected_positions.pop()

# Print the selected positions
print(selected_positions)