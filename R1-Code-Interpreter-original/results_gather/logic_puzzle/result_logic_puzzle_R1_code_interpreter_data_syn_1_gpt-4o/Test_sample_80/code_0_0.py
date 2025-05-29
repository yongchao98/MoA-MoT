# Define the grid
grid = [
    [1, -5, -7, 1, -2, 1],
    [-9, 3, -4, -3, 9, 6],
    [2, 12, -1, -9, 5, 11],
    [10, 11, 11, 10, 2, 12],
    [2, 9, 4, 6, 3, 7],
    [1, -5, 5, 6, 1, 1]
]

# Initialize the selected positions with the given starting point
selected_positions = [(0, 5)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Try to select 12 numbers
for r in range(6):
    for c in range(6):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 12 and is_valid_selection(selected_positions):
                break
            if len(selected_positions) > 12 or not is_valid_selection(selected_positions):
                selected_positions.pop()
    if len(selected_positions) == 12:
        break

# Print the selected positions
print(selected_positions)