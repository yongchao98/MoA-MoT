import itertools

# Define the grid
grid = [
    [2, 1, 2],
    [3, 1, 3],
    [4, 2, 2]
]

# The first selected number is at position (0, 2)
selected_positions = [(0, 2)]
selected_numbers = {grid[0][2]}

# Function to check if the selection is valid
def is_valid_selection(positions):
    row_sums = [0] * 3
    col_sums = [0] * 3
    numbers = set()
    
    for r, c in positions:
        number = grid[r][c]
        if number in numbers:
            return False
        numbers.add(number)
        row_sums[r] += number
        col_sums[c] += number
    
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try all combinations of 3 additional positions
for additional_positions in itertools.combinations([(r, c) for r in range(3) for c in range(3) if (r, c) != (0, 2)], 3):
    positions = selected_positions + list(additional_positions)
    if is_valid_selection(positions):
        print(positions)
        break