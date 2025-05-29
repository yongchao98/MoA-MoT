grid = [
    [-10, -10, 11, 3, -11, -12],
    [5, 2, 8, -1, 1, -4],
    [-12, 11, 11, 10, 8, 10],
    [7, -8, -9, -5, 12, 2],
    [9, 12, 5, 7, 11, 2],
    [-8, 9, -6, 1, 5, -9]
]

# Initial selection
selected_positions = [(0, 5)]
selected_numbers = [-12]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_count = [0] * 6
    col_count = [0] * 6
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selected_positions:
        row_count[r] += 1
        col_count[c] += 1
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for i in range(6):
        if row_count[i] > 0 and row_neg_count[i] % 2 != 0:
            return False
        if col_count[i] > 0 and col_neg_count[i] % 2 != 0:
            return False
    
    return True

# Try to select 12 numbers
for r in range(6):
    for c in range(6):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            selected_numbers.append(grid[r][c])
            if len(selected_positions) == 12 and is_valid_selection(selected_positions):
                break
            if len(selected_positions) > 12 or not is_valid_selection(selected_positions):
                selected_positions.pop()
                selected_numbers.pop()
    if len(selected_positions) == 12:
        break

print(selected_positions)