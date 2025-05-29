# Define the grid
grid = [
    [-13, 6, 12, 14, 12, -3, -7],
    [13, -13, 1, 13, 13, 2, -10],
    [2, 6, -9, 13, 1, 10, 12],
    [8, -6, 13, 8, 14, 9, 13],
    [3, 10, -14, 10, -6, 5, 11],
    [-9, 11, 4, 5, 7, -11, -10],
    [-5, 1, 4, 13, 14, 2, 12]
]

# Start with the given number -7 at position (0,6)
selected_positions = [(0, 6)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 7
    col_neg_count = [0] * 7
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    # Check if all rows and columns have even negative counts
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Backtracking function to find a valid selection
def find_selection(selected_positions):
    if len(selected_positions) == 16:
        return selected_positions
    for r in range(7):
        for c in range(7):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if is_valid_selection(selected_positions):
                    result = find_selection(selected_positions)
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find a valid selection
result = find_selection(selected_positions)

# Print the selected positions
print(result)