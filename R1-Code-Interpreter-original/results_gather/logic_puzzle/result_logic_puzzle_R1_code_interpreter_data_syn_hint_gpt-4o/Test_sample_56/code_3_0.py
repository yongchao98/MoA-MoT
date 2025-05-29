# Given grid
grid = [
    [-5, 7, -7, 6, -1],
    [-2, 2, 8, 6, -3],
    [3, 4, 4, 3, 8],
    [-8, -2, 5, 8, 2],
    [-2, 6, 8, 2, -1]
]

# Start with the given number -1 at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    return all(count % 2 == 0 for count in row_count) and all(count % 2 == 0 for count in col_count)

# Systematic selection of numbers
def find_valid_selection():
    # Pre-select numbers to ensure even negatives in each row and column
    # Manually select numbers to satisfy the constraints
    positions = [
        (0, 4),  # Given
        (0, 0), (0, 1), (0, 3),  # Row 0
        (1, 0), (1, 1), (1, 2),  # Row 1
        (2, 0), (2, 1),  # Row 2
    ]
    if is_valid_selection(positions):
        return positions
    return None

# Find a valid selection
result = find_valid_selection()
print(result)