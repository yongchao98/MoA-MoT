import itertools

# Given grid
grid = [
    [-5, 7, -7, 6, -1],
    [-2, 2, 8, 6, -3],
    [3, 4, 4, 3, 8],
    [-8, -2, 5, 8, 2],
    [-2, 6, 8, 2, -1]
]

# Start with the given number -1 at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    return all(count % 2 == 0 for count in row_count) and all(count % 2 == 0 for count in col_count)

# Try to find a valid selection of 9 numbers
for positions in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    positions = selected_positions + list(positions)
    if is_valid_selection(positions):
        print(positions)
        break