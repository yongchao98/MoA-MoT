# Define the grid
grid = [
    [6, 1, 5, -9, 3, -10],
    [-9, -4, -9, 2, 7, 11],
    [9, -8, 6, 12, 12, -2],
    [9, -2, 11, -6, 8, 12],
    [4, 1, 2, -2, -5, 10],
    [3, -5, 6, -1, 9, 12]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    # Check if each row and column has an even number of negative numbers
    for i in range(6):
        if row_neg_count[i] % 2 != 0 or col_neg_count[i] % 2 != 0:
            return False
    return True

# Try to select 12 numbers
def select_numbers():
    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if len(selected_positions) == 12 and is_valid_selection(selected_positions):
                    return selected_positions
                if len(selected_positions) < 12:
                    result = select_numbers()
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find the selection
result = select_numbers()
print(result)