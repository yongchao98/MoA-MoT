# Define the grid
grid = [
    [4, 4, 1],
    [3, 2, 1],
    [2, 2, 4]
]

# Start with the given selected number
selected_positions = [(0, 2)]
selected_numbers = [grid[0][2]]

# Function to check if adding a number at a position is valid
def is_valid_selection(r, c, selected_positions, selected_numbers):
    # Check if the number is already selected
    if grid[r][c] in selected_numbers:
        return False
    
    # Calculate the sum for the row and column
    row_sum = sum(grid[r][j] for i, j in selected_positions if i == r) + grid[r][c]
    col_sum = sum(grid[i][c] for i, j in selected_positions if j == c) + grid[r][c]
    
    # Check if the sum is less than or equal to 4
    return row_sum <= 4 and col_sum <= 4

# Try to select 3 more numbers
for r in range(3):
    for c in range(3):
        if (r, c) not in selected_positions and is_valid_selection(r, c, selected_positions, selected_numbers):
            selected_positions.append((r, c))
            selected_numbers.append(grid[r][c])
            if len(selected_positions) == 4:
                break
    if len(selected_positions) == 4:
        break

# Print the selected positions
print(selected_positions)