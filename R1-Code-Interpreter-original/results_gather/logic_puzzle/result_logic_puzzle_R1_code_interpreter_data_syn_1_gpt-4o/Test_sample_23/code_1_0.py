import itertools

# Define the grid
grid = [
    [3, 4, 2],
    [1, 3, 4],
    [4, 4, 3]
]

# Pre-selected number
selected_positions = [(0, 2)]
selected_numbers = [grid[0][2]]

# All possible positions except the pre-selected one
all_positions = [(i, j) for i in range(3) for j in range(3) if (i, j) != (0, 2)]

# Function to check if the selection is valid
def is_valid_selection(positions):
    row_sums = [0] * 3
    col_sums = [0] * 3
    numbers = set()
    
    for r, c in positions:
        num = grid[r][c]
        if num in numbers:
            return False
        numbers.add(num)
        row_sums[r] += num
        col_sums[c] += num
    
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try all combinations of 3 positions from the remaining positions
for combination in itertools.combinations(all_positions, 3):
    positions = selected_positions + list(combination)
    if is_valid_selection(positions):
        print(positions)
        break