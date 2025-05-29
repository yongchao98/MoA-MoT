from itertools import combinations

# Define the grid
grid = [
    [-5, 7, -7, 6, -1],
    [-2, 2, 8, 6, -3],
    [3, 4, 4, 3, 8],
    [-8, -2, 5, 8, 2],
    [-2, 6, 8, 2, -1]
]

# Pre-select the given number
selected_positions = [(0, 4)]

# Function to check if a selection is valid
def is_valid_selection(positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in positions:
        row_count[r] += 1
        col_count[c] += 1
    
    # Check if each row and column has an even number of negative numbers
    for r in range(5):
        if row_count[r] % 2 != 0:
            return False
    for c in range(5):
        if col_count[c] % 2 != 0:
            return False
    return True

# Try to find a valid combination
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)]
for combo in combinations(all_positions, 8):
    current_selection = selected_positions + list(combo)
    if is_valid_selection(current_selection):
        print(current_selection)
        break