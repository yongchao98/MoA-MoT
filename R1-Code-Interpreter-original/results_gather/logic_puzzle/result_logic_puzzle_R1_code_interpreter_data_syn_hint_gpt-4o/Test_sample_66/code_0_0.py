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

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    row_neg_count = [0] * 5
    col_neg_count = [0] * 5
    
    for r, c in selected_positions:
        row_count[r] += 1
        col_count[c] += 1
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for i in range(5):
        if row_count[i] > 0 and row_neg_count[i] % 2 != 0:
            return False
        if col_count[i] > 0 and col_neg_count[i] % 2 != 0:
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
            elif len(selected_positions) > 10:
                break

# Print the selected positions
print(selected_positions)