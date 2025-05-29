import itertools

# Define the grid
grid = [
    [3, 4, 1],
    [4, 4, 2],
    [3, 3, 1]
]

# Pre-select the number at (0,2)
selected_positions = [(0, 2)]
selected_numbers = [grid[0][2]]

# All possible positions except the pre-selected one
all_positions = [(i, j) for i in range(3) for j in range(3) if (i, j) != (0, 2)]

# Function to check if a selection is valid
def is_valid_selection(positions):
    # Check if all numbers are different
    numbers = [grid[i][j] for i, j in positions]
    if len(set(numbers)) != len(numbers):
        return False
    
    # Check row and column sums
    row_sums = [0] * 3
    col_sums = [0] * 3
    for i, j in positions:
        row_sums[i] += grid[i][j]
        col_sums[j] += grid[i][j]
    
    if all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums):
        return True
    return False

# Try all combinations of 3 positions from the remaining positions
for positions in itertools.combinations(all_positions, 3):
    full_selection = selected_positions + list(positions)
    if is_valid_selection(full_selection):
        # Convert positions to the required format
        result = ', '.join(f"{i} {j}" for i, j in full_selection)
        print(result)
        break