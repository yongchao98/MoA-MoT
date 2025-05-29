# The grid
grid = [
    [5, -2, 2, -8, -6],
    [1, -2, 7, -1, 8],
    [4, -7, 6, -7, -3],
    [8, -7, 6, 8, 8],
    [-4, -5, 4, 6, 3]
]

# Pre-selected number
pre_selected = (0, 4)

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    # Check rows
    for i in range(5):
        row_product = 1
        for j in range(5):
            if (i, j) in selected_positions:
                row_product *= grid[i][j]
        if row_product <= 0:
            return False
    
    # Check columns
    for j in range(5):
        col_product = 1
        for i in range(5):
            if (i, j) in selected_positions:
                col_product *= grid[i][j]
        if col_product <= 0:
            return False
    
    return True

# Generate all possible combinations of 9 positions including the pre-selected one
all_positions = [(i, j) for i in range(5) for j in range(5)]
all_positions.remove(pre_selected)

# Try all combinations
for combination in itertools.combinations(all_positions, 8):
    selected_positions = set(combination) | {pre_selected}
    if is_valid_selection(selected_positions):
        # Convert to the required format
        result = ', '.join(f"{i} {j}" for i, j in selected_positions)
        print(f"<<<{result}>>>")
        break