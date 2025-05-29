# Define the grid
grid = [
    [4, 7, -7, 6, 2],
    [-6, -3, 3, 10, 2],
    [9, 4, 10, 5, 4],
    [-4, 4, 8, 5, 7],
    [8, 1, 10, -2, 10]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 4)]

# Function to check if adding a number keeps the product positive
def is_valid_selection(selected_positions, new_position):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    r, c = new_position
    if grid[r][c] < 0:
        row_count[r] += 1
        col_count[c] += 1
    # Check if the number of negative numbers in each row and column is even
    for i in range(5):
        if row_count[i] % 2 != 0 or col_count[i] % 2 != 0:
            return False
    return True

# Select numbers
for r in range(5):
    for c in range(5):
        if (r, c) != (0, 4) and len(selected_positions) < 10:
            if is_valid_selection(selected_positions, (r, c)):
                selected_positions.append((r, c))

# Print the selected positions
print(selected_positions)