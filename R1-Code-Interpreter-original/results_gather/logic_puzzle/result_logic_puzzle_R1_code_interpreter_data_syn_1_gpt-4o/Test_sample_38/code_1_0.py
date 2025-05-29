import itertools

# Define the grid
grid = [
    [3, -4, 3, -1],
    [8, -4, -2, 1],
    [7, 4, -4, 2],
    [5, 4, 5, 3]
]

# Pre-select the given number
selected_positions = [(0, 3)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_count = [0] * 4
    col_count = [0] * 4
    row_neg_count = [0] * 4
    col_neg_count = [0] * 4
    
    for r, c in selected_positions:
        row_count[r] += 1
        col_count[c] += 1
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for i in range(4):
        if row_count[i] > 0 and row_neg_count[i] % 2 != 0:
            return False
        if col_count[i] > 0 and col_neg_count[i] % 2 != 0:
            return False
    
    return True

# Try all combinations of 5 more positions
all_positions = [(r, c) for r in range(4) for c in range(4) if (r, c) != (0, 3)]
for combination in itertools.combinations(all_positions, 5):
    selected = selected_positions + list(combination)
    if is_valid_selection(selected):
        print(selected)
        break