import itertools

# Define the grid
grid = [
    [3, 1, -1, -6, 7],
    [7, -4, -3, 1, 6],
    [1, 3, -2, 5, 3],
    [6, 8, -7, -8, 5],
    [1, 4, 5, 1, -7]
]

# Pre-select the given number
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
        row_product = 1
        for c in range(5):
            if (r, c) in positions:
                row_product *= grid[r][c]
        if row_product <= 0:
            return False
    
    for c in range(5):
        col_product = 1
        for r in range(5):
            if (r, c) in positions:
                col_product *= grid[r][c]
        if col_product <= 0:
            return False
    
    return True

# Try all combinations of 8 additional numbers
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)]
for combination in itertools.combinations(all_positions, 8):
    current_selection = selected_positions + list(combination)
    if is_valid_selection(current_selection):
        print(current_selection)
        break