# Define the grid
grid = [
    [-7, 8, 9, 3, 2, -9],
    [-5, -8, 6, 10, 12, 2],
    [-4, 5, 12, 7, 6, 12],
    [1, 7, 7, -2, -5, 6],
    [-10, 6, -5, 9, 2, -8],
    [10, -10, 2, 1, -1, -8]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if a selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    # Check if all rows and columns have even number of negatives
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