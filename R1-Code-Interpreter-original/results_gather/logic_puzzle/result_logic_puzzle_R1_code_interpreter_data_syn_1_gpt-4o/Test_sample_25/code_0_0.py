# Define the grid
grid = [
    [2, 1, 1],
    [4, 2, 3],
    [3, 3, 3]
]

# Start with the given selected number
selected_positions = [(0, 2)]
selected_numbers = [grid[0][2]]

# Function to check if adding a number keeps the sum in row and column <= 4
def is_valid_selection(selected_positions, new_position, grid):
    row, col = new_position
    # Calculate row and column sums
    row_sum = sum(grid[r][col] for r, c in selected_positions if r == row) + grid[row][col]
    col_sum = sum(grid[r][c] for r, c in selected_positions if c == col) + grid[row][col]
    return row_sum <= 4 and col_sum <= 4

# Try to select 3 more numbers
for r in range(3):
    for c in range(3):
        if (r, c) not in selected_positions and grid[r][c] not in selected_numbers:
            if is_valid_selection(selected_positions, (r, c), grid):
                selected_positions.append((r, c))
                selected_numbers.append(grid[r][c])
            if len(selected_positions) == 4:
                break
    if len(selected_positions) == 4:
        break

# Print the selected positions
print(selected_positions)