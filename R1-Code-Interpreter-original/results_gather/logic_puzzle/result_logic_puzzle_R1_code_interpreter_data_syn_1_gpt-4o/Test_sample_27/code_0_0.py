import itertools

# Define the grid
grid = [
    [3, 4, 3],
    [3, 1, 2],
    [1, 2, 3]
]

# Pre-selected position
pre_selected = (0, 2)

# Function to check if a selection is valid
def is_valid_selection(selection):
    # Check if all numbers are different
    selected_numbers = [grid[r][c] for r, c in selection]
    if len(set(selected_numbers)) != len(selected_numbers):
        return False
    
    # Check row and column sums
    row_sums = [0] * 3
    col_sums = [0] * 3
    for r, c in selection:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    if any(s > 4 for s in row_sums) or any(s > 4 for s in col_sums):
        return False
    
    return True

# Generate all possible selections of 4 positions including the pre-selected one
all_positions = [(r, c) for r in range(3) for c in range(3)]
all_positions.remove(pre_selected)

# Try all combinations
for selection in itertools.combinations(all_positions, 3):
    selection = [pre_selected] + list(selection)
    if is_valid_selection(selection):
        # Format the output as required
        result = ', '.join(f"{r} {c}" for r, c in selection)
        print(result)
        break