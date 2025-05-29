import itertools

# Define the grid
grid = [
    [3, 1, 3],
    [1, 1, 2],
    [4, 3, 4]
]

# Pre-select the number at (0, 2)
selected_positions = [(0, 2)]
selected_numbers = {grid[0][2]}

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

# Try to find a valid selection
for positions in itertools.combinations([(r, c) for r in range(3) for c in range(3) if (r, c) != (0, 2)], 3):
    positions = selected_positions + list(positions)
    if is_valid_selection(positions):
        print(positions)
        break