# Define the grid
grid = [
    [-9, 8, 3, 13, 5, -1, 6],
    [13, -6, -4, 6, 13, 6, -10],
    [7, -6, -2, 7, 6, -9, 3],
    [6, -5, 1, 4, -10, 2, 9],
    [-9, 14, 12, 1, -10, 9, -2],
    [-3, -2, 14, 14, 2, -8, 5],
    [-2, -1, -8, -3, 7, -6, 9]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 6)]

# Function to check if adding a number keeps the product positive
def is_valid_selection(selected_positions, new_position):
    row_counts = [0] * 7
    col_counts = [0] * 7
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_counts[r] += 1
            col_counts[c] += 1
    r, c = new_position
    if grid[r][c] < 0:
        row_counts[r] += 1
        col_counts[c] += 1
    # Check if all rows and columns have even negative counts or are all positive
    for count in row_counts + col_counts:
        if count % 2 != 0:
            return False
    return True

# Select numbers
for r in range(7):
    for c in range(7):
        if (r, c) != (0, 6) and len(selected_positions) < 16:
            if is_valid_selection(selected_positions, (r, c)):
                selected_positions.append((r, c))

# Format the output
output = ', '.join(f"{r} {c}" for r, c in selected_positions)
print(output)