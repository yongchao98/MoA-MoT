import itertools

# Define the grid
grid = [
    [3, 3, 1],
    [4, 3, 2],
    [4, 4, 3]
]

# Pre-selected position
selected_positions = [(0, 2)]

# Function to check if the sum constraints are satisfied
def is_valid_selection(positions):
    row_sums = [0, 0, 0]
    col_sums = [0, 0, 0]
    
    for r, c in positions:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Generate all possible combinations of 3 additional positions
all_positions = [(r, c) for r in range(3) for c in range(3) if (r, c) != (0, 2)]
combinations = itertools.combinations(all_positions, 3)

# Check each combination
for combo in combinations:
    positions = selected_positions + list(combo)
    if is_valid_selection(positions):
        print(positions)
        break