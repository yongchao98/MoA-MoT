import itertools

# Define the grid
grid = [
    [3, 1, 1],
    [4, 1, 4],
    [4, 2, 4]
]

# Pre-select the number at (0, 2)
selected_positions = [(0, 2)]
selected_numbers = {grid[0][2]}

# Function to check if the current selection is valid
def is_valid_selection(positions):
    row_sums = [0] * 3
    col_sums = [0] * 3
    for r, c in positions:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try to find the remaining 3 positions
for positions in itertools.combinations([(r, c) for r in range(3) for c in range(3) if (r, c) != (0, 2)], 3):
    # Check if all numbers are different
    numbers = {grid[r][c] for r, c in positions}
    if len(numbers) == 3 and not numbers.intersection(selected_numbers):
        # Check if the selection is valid
        if is_valid_selection(selected_positions + list(positions)):
            selected_positions += list(positions)
            break

# Output the selected positions
print(selected_positions)