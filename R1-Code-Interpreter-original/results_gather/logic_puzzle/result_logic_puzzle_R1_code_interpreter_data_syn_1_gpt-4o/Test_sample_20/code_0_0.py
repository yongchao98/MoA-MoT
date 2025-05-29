import itertools

# Define the grid
grid = [
    [3, 1, 2],
    [2, 4, 4],
    [3, 2, 1]
]

# Pre-selected number
selected_positions = [(0, 2)]
selected_numbers = [grid[0][2]]

# All possible positions except the pre-selected one
all_positions = [(i, j) for i in range(3) for j in range(3) if (i, j) != (0, 2)]

# Function to check if the selection is valid
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
    
    if any(s > 4 for s in row_sums) or any(s > 4 for s in col_sums):
        return False
    
    return True

# Try all combinations of 3 positions from the remaining positions
for combination in itertools.combinations(all_positions, 3):
    # Add the pre-selected position
    current_selection = selected_positions + list(combination)
    if is_valid_selection(current_selection):
        # If valid, print the positions
        print(current_selection)
        break