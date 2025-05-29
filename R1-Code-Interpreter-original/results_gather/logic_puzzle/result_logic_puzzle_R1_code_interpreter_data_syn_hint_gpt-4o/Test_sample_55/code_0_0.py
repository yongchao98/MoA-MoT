import itertools

# Define the grid
grid = [
    [7, -4, -6, 3, -3],
    [4, 8, 5, 4, 3],
    [2, -4, 7, -8, -5],
    [6, -8, -1, 4, -8],
    [8, 8, -3, -7, -2]
]

# Start with the given number at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in positions:
        row_count[r] += 1
        col_count[c] += 1
    
    # Check if each row and column has an even number of negative numbers
    for r in range(5):
        row_neg_count = sum(1 for c in range(5) if (r, c) in positions and grid[r][c] < 0)
        if row_count[r] > 0 and row_neg_count % 2 != 0:
            return False
    
    for c in range(5):
        col_neg_count = sum(1 for r in range(5) if (r, c) in positions and grid[r][c] < 0)
        if col_count[c] > 0 and col_neg_count % 2 != 0:
            return False
    
    return True

# Try to find a valid selection of 9 numbers
for positions in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    positions = selected_positions + list(positions)
    if is_valid_selection(positions):
        print(positions)
        break