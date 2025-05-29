import itertools

# Define the grid
grid = [
    [-2, 3, 7, 8, 8],
    [-5, 1, 5, 8, 6],
    [7, -3, -1, -4, 2],
    [5, -4, -7, 3, 3],
    [2, -4, 3, 1, 6]
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

# Generate all possible combinations of 8 additional positions
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != initial_position]
for combination in itertools.combinations(all_positions, 8):
    selected_positions = [initial_position] + list(combination)
    if is_valid_selection(selected_positions):
        print(selected_positions)
        break