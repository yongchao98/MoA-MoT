from itertools import combinations

# Define the grid
grid = [
    [1, 3, 2],
    [4, 1, 3],
    [3, 3, 4]
]

# Pre-select the number 2 at (0,2)
selected_positions = [(0, 2)]
selected_numbers = {2}

# Function to check if a selection is valid
def is_valid_selection(positions):
    row_sums = [0, 0, 0]
    col_sums = [0, 0, 0]
    numbers = set(selected_numbers)
    
    for r, c in positions:
        num = grid[r][c]
        if num in numbers:
            return False
        numbers.add(num)
        row_sums[r] += num
        col_sums[c] += num
    
    # Add the pre-selected number to the sums
    row_sums[0] += 2
    col_sums[2] += 2
    
    if all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums):
        return True
    return False

# Generate all possible combinations of 3 positions from the grid
all_positions = [(r, c) for r in range(3) for c in range(3) if (r, c) != (0, 2)]
result = None
for combo in combinations(all_positions, 3):
    if is_valid_selection(combo):
        result = selected_positions + list(combo)
        break

# Print the result if found
if result:
    print(result)
else:
    print("No valid selection found")