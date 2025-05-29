import itertools

# Define the grid
grid = [
    [-1, -6, 5, 5, -8],
    [1, -5, 5, -1, -4],
    [5, 8, 6, 5, 3],
    [7, 5, 2, 3, 5],
    [6, 4, 6, 7, 3]
]

# Initial selection
selected_positions = [(0, 4)]

# Function to check if the selection is valid
def is_valid_selection(selected_positions):
    # Check row and column products
    row_counts = [0] * 5
    col_counts = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_counts[r] += 1
            col_counts[c] += 1
    # Check if all row and column counts are even
    return all(count % 2 == 0 for count in row_counts) and all(count % 2 == 0 for count in col_counts)

# Try to find a valid selection
for positions in itertools.combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 8):
    selected_positions = [(0, 4)] + list(positions)
    if is_valid_selection(selected_positions):
        break

# Print the selected positions
print(selected_positions)