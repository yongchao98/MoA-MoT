import itertools

# Define the grid
grid = [
    [7, -4, -6, 3, -3],
    [4, 8, 5, 4, 3],
    [2, -4, 7, -8, -5],
    [6, -8, -1, 4, -8],
    [8, 8, -3, -7, -2]
]

# Pre-select the given number
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
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
for positions in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    selected_positions = [(0, 4)] + list(positions)
    if is_valid_selection(selected_positions):
        print(selected_positions)
        break