# Define the grid
grid = [
    [2, 1],
    [1, 3]
]

# Already selected number
selected_positions = [(0, 1)]

# Function to check if the selection is valid
def is_valid_selection(selected_positions):
    row_sums = [0, 0]
    col_sums = [0, 0]
    
    for r, c in selected_positions:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try all combinations of selecting two more numbers
from itertools import combinations

# Remaining positions to select from
remaining_positions = [(0, 0), (1, 0), (1, 1)]

# Check all combinations of selecting two more numbers
for combo in combinations(remaining_positions, 2):
    current_selection = selected_positions + list(combo)
    if is_valid_selection(current_selection):
        print(current_selection)
        break