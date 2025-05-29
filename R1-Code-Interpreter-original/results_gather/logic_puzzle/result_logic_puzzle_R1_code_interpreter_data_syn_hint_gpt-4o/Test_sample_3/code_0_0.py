# Define the grid
grid = [
    [2, 1],
    [3, 3]
]

# Pre-selected number
selected = [(0, 1)]

# Function to check if the selection is valid
def is_valid_selection(selection):
    row_sums = [0, 0]
    col_sums = [0, 0]
    
    for r, c in selection:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try all combinations of selecting two more numbers
from itertools import combinations

# All possible positions except the pre-selected one
positions = [(0, 0), (1, 0), (1, 1)]

# Check all combinations
for comb in combinations(positions, 2):
    selection = selected + list(comb)
    if is_valid_selection(selection):
        print(selection)
        break
else:
    print("No valid selection found")