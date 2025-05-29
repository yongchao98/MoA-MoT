import itertools

# The grid
grid = [
    [1, 4, 4],
    [2, 1, 4],
    [2, 3, 3]
]

# The first selected number is 4 at (0,2)
selected_positions = [(0, 2)]
selected_numbers = {4}

# Possible positions to select from
possible_positions = [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]

# Try all combinations of 3 positions from the possible positions
for positions in itertools.combinations(possible_positions, 3):
    numbers = {grid[r][c] for r, c in positions}
    if len(numbers) == 3 and not numbers.intersection(selected_numbers):
        # Check row and column sums
        row_sums = [0, 0, 0]
        col_sums = [0, 0, 0]
        
        # Add the first selected number
        row_sums[0] += 4
        col_sums[2] += 4
        
        # Add the new selected numbers
        for r, c in positions:
            row_sums[r] += grid[r][c]
            col_sums[c] += grid[r][c]
        
        # Check if all sums are <= 4
        if all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums):
            # Found a valid selection
            selected_positions.extend(positions)
            break

# Output the selected positions
print(selected_positions)