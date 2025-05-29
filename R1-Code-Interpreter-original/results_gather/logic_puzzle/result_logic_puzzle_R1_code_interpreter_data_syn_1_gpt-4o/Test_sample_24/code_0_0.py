import itertools

# Define the grid
grid = [
    [1, 3, 4],
    [3, 3, 3],
    [2, 2, 2]
]

# Pre-selected position
selected_positions = [(0, 2)]

# Possible positions to select from
possible_positions = [(i, j) for i in range(3) for j in range(3) if (i, j) != (0, 2)]

# Function to check if a selection is valid
def is_valid_selection(positions):
    row_sums = [0] * 3
    col_sums = [0] * 3
    selected_numbers = set()
    
    for pos in positions:
        row, col = pos
        num = grid[row][col]
        if num in selected_numbers:
            return False
        selected_numbers.add(num)
        row_sums[row] += num
        col_sums[col] += num
    
    # Check if all row and column sums are <= 4
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try all combinations of 3 positions from possible_positions
for combination in itertools.combinations(possible_positions, 3):
    if is_valid_selection(selected_positions + list(combination)):
        result_positions = selected_positions + list(combination)
        break

# Format the result
result = ', '.join(f"{r} {c}" for r, c in result_positions)
print(result)