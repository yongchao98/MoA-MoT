from itertools import combinations

# Grid representation
grid = [
    [2, 3],
    [3, 1]
]

# Pre-selected number at (0,1)
selected_positions = [(0, 1)]

# Possible positions to select from
possible_positions = [(0, 0), (1, 0), (1, 1)]

# Check all combinations of two additional selections
for combo in combinations(possible_positions, 2):
    # Combine with the pre-selected position
    selected = selected_positions + list(combo)
    
    # Calculate sums for rows and columns
    row_sums = [0, 0]
    col_sums = [0, 0]
    
    for r, c in selected:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    # Check if all sums are <= 4
    if all(s <= 4 for s in row_sums + col_sums):
        print(selected)
        break