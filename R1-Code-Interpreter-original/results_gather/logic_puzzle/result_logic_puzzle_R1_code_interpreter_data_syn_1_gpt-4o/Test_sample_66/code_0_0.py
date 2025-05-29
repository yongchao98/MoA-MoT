# Define the grid
grid = [
    [6, -8, -3, 9, 8],
    [4, 8, 9, 3, -8],
    [9, 4, 5, 2, 1],
    [-4, -2, -4, 2, 7],
    [1, -3, 9, 1, -10]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    # Check if all rows and columns have an even number of negative numbers
    return all(count % 2 == 0 for count in row_count) and all(count % 2 == 0 for count in col_count)

# Try to select 9 more numbers
for r in range(5):
    for c in range(5):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 10 and is_product_positive(selected_positions, grid):
                break
            if len(selected_positions) > 10 or not is_product_positive(selected_positions, grid):
                selected_positions.pop()
    if len(selected_positions) == 10:
        break

# Print the selected positions
print(selected_positions)