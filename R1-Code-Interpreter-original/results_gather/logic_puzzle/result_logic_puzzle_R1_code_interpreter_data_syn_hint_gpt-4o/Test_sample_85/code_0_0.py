# Define the grid
grid = [
    [9, -11, -5, 1, 7, 9],
    [-12, 10, 3, -5, 3, 6],
    [12, -12, -5, 4, 11, 11],
    [4, 12, 8, -6, 7, 7],
    [12, 3, 4, 7, 12, 4],
    [5, 12, 6, 1, 2, 9]
]

# Start with the given number
selected_positions = [(0, 5)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    # Check if all rows and columns have even negative counts
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Try to select 12 numbers
for r in range(6):
    for c in range(6):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            if len(selected_positions) == 12:
                if is_valid_selection(selected_positions):
                    print(selected_positions)
                    break
            else:
                if not is_valid_selection(selected_positions):
                    selected_positions.pop()
    if len(selected_positions) == 12:
        break