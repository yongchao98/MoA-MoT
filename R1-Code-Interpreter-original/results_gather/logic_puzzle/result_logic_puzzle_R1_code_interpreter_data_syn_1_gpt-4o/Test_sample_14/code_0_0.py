# Define the grid
grid = [
    [2, 1],
    [3, 1]
]

# Pre-selected position
selected_positions = [(0, 1)]

# Possible positions to select from
possible_positions = [(0, 0), (1, 0), (1, 1)]

# Function to check if the selection is valid
def is_valid_selection(selected_positions):
    row_sums = [0, 0]
    col_sums = [0, 0]
    
    for r, c in selected_positions:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    return all(s <= 4 for s in row_sums + col_sums)

# Try all combinations of two additional selections
from itertools import combinations

for additional_selection in combinations(possible_positions, 2):
    current_selection = selected_positions + list(additional_selection)
    if is_valid_selection(current_selection):
        print(current_selection)
        break