# Define the grid
grid = [
    [3, 1],
    [2, 3]
]

# Pre-selected number
selected_positions = [(0, 1)]  # (0,1) is already selected

# Possible positions to select from
possible_positions = [(0, 0), (1, 0), (1, 1)]

# Check all combinations of selecting two more numbers
from itertools import combinations

# Function to check if the selection is valid
def is_valid_selection(selected_positions):
    # Calculate sums for rows and columns
    row_sums = [0, 0]
    col_sums = [0, 0]
    
    for r, c in selected_positions:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    # Check if all sums are <= 4
    return all(s <= 4 for s in row_sums + col_sums)

# Try all combinations of two more selections
for pos1, pos2 in combinations(possible_positions, 2):
    current_selection = selected_positions + [pos1, pos2]
    if is_valid_selection(current_selection):
        # Convert to the required format
        result = ', '.join(f"{r} {c}" for r, c in current_selection)
        print(f"<<<{result}>>>")
        break