import itertools

# The given grid
grid = [
    [8, -3, 5, 1, 3],
    [8, 5, -6, 3, 4],
    [5, -8, 7, -1, -8],
    [-2, -3, -4, -7, 1],
    [-1, 8, -8, -4, 8]
]

# Initial selected position
initial_position = (0, 4)

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    row_neg_count = [0] * 5
    col_neg_count = [0] * 5
    
    for r, c in selected_positions:
        row_count[r] += 1
        col_count[c] += 1
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for i in range(5):
        if row_count[i] > 0 and (row_neg_count[i] % 2 != 0):
            return False
        if col_count[i] > 0 and (col_neg_count[i] % 2 != 0):
            return False
    
    return True

# Generate all possible combinations of 9 positions including the initial one
all_positions = [(r, c) for r in range(5) for c in range(5)]
all_positions.remove(initial_position)

for combination in itertools.combinations(all_positions, 8):
    selected_positions = [initial_position] + list(combination)
    if is_valid_selection(selected_positions):
        print(selected_positions)
        break