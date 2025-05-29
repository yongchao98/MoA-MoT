# Define the grid
grid = [
    [2, 8, -12, -1, 9, 7],
    [2, 7, 6, 3, 1, -1],
    [9, 12, 9, 8, -10, 2],
    [6, 10, -5, -12, 5, -7],
    [-12, -8, 3, 8, 12, 4],
    [-5, -9, 11, -12, 5, -8]
]

# Initialize selected positions with the given number
selected_positions = [(0, 5)]
selected_numbers = [7]

# Function to check if adding a number maintains positive product
def is_valid_selection(selected_positions, row, col):
    row_neg_count = sum(1 for r, c in selected_positions if r == row and grid[r][c] < 0)
    col_neg_count = sum(1 for r, c in selected_positions if c == col and grid[r][c] < 0)
    current_number = grid[row][col]
    if current_number < 0:
        row_neg_count += 1
        col_neg_count += 1
    return row_neg_count % 2 == 0 and col_neg_count % 2 == 0

# Select numbers
for r in range(6):
    for c in range(6):
        if (r, c) not in selected_positions:
            if is_valid_selection(selected_positions, r, c):
                selected_positions.append((r, c))
                selected_numbers.append(grid[r][c])
                if len(selected_positions) == 12:
                    break
    if len(selected_positions) == 12:
        break

# Print the selected positions
print(selected_positions)