import itertools

# The given grid
grid = [
    [5, 1, 10, 11, 1, 1],
    [7, 9, -1, 12, -9, 8],
    [3, 10, 1, -1, -11, 1],
    [3, -10, -4, -9, 3, 6],
    [3, -2, 3, -10, 2, 3],
    [-9, -3, 8, 4, 6, 6]
]

# Initial selected position
initial_position = (0, 5)

# Function to check if a selection is valid
def is_valid_selection(selection):
    # Check row and column products
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selection:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    # Check if all rows and columns have even negative counts
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Generate all possible selections of 12 numbers including the initial position
all_positions = [(r, c) for r in range(6) for c in range(6)]
all_positions.remove(initial_position)

# Try all combinations
for selection in itertools.combinations(all_positions, 11):
    full_selection = [initial_position] + list(selection)
    if is_valid_selection(full_selection):
        # Convert to the required format
        result = ', '.join(f"{r} {c}" for r, c in full_selection)
        print(result)
        break