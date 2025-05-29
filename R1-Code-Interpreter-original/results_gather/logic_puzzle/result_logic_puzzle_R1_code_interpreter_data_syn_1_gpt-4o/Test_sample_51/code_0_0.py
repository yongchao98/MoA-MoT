import itertools

# Define the grid
grid = [
    [3, -7, 5, 6, 6],
    [1, 8, -2, 6, 7],
    [-6, 1, -7, 7, -6],
    [4, -1, -8, -5, -3],
    [-3, 2, -8, -6, 5]
]

# Initial selected position
initial_position = (0, 4)

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    # Check if each row and column has an even number of negative numbers
    for count in row_count + col_count:
        if count % 2 != 0:
            return False
    return True

# Generate all possible combinations of 9 positions including the initial position
all_positions = [(r, c) for r in range(5) for c in range(5)]
all_positions.remove(initial_position)

# Iterate over combinations
for combination in itertools.combinations(all_positions, 8):
    selected_positions = [initial_position] + list(combination)
    if is_valid_selection(selected_positions):
        # Convert positions to the required format
        result = ', '.join(f"{r} {c}" for r, c in selected_positions)
        print(result)
        break