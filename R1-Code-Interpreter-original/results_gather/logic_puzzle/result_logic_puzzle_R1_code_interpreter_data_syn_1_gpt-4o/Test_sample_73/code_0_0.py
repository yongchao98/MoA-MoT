import itertools

# The grid
grid = [
    [-8, 2, -9, 9, -6],
    [4, 2, -2, 7, 7],
    [-7, 6, 8, -4, -5],
    [10, -1, -7, 4, 1],
    [-8, -9, 3, 7, 10]
]

# Initial selection
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    # Check if each row and column has an even number of negative numbers
    for count in row_count + col_count:
        if count % 2 != 0:
            return False
    return True

# Try to find a valid selection
for positions in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 9):
    selected_positions = [(0, 4)] + list(positions)
    if is_valid_selection(selected_positions):
        break

# Output the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))