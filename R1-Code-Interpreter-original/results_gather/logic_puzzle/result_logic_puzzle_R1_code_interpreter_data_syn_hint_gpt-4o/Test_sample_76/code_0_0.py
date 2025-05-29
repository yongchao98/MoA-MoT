# Define the grid
grid = [
    [-8, 10, -8, -3, 2],
    [4, 10, 5, 2, 9],
    [-5, -4, 5, 9, -2],
    [1, 1, 5, -4, -9],
    [-9, -5, 1, 10, 3]
]

# Start with the given number at (0,4)
selected_positions = [(0, 4)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    row_product = [1] * 5
    col_product = [1] * 5
    
    for r, c in selected_positions:
        row_count[r] += 1
        col_count[c] += 1
        row_product[r] *= grid[r][c]
        col_product[c] *= grid[r][c]
    
    for i in range(5):
        if row_count[i] > 0 and row_product[i] <= 0:
            return False
        if col_count[i] > 0 and col_product[i] <= 0:
            return False
    
    return True

# Try to select 10 numbers
for r in range(5):
    for c in range(5):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 10:
                if is_valid_selection(selected_positions):
                    break
                else:
                    selected_positions.pop()
            elif not is_valid_selection(selected_positions):
                selected_positions.pop()
    if len(selected_positions) == 10:
        break

# Print the selected positions
print(selected_positions)