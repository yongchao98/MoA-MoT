import itertools

# Define the grid
grid = [
    [-10, -10, 11, 3, -11, -12],
    [5, 2, 8, -1, 1, -4],
    [-12, 11, 11, 10, 8, 10],
    [7, -8, -9, -5, 12, 2],
    [9, 12, 5, 7, 11, 2],
    [-8, 9, -6, 1, 5, -9]
]

# Start with the given position
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_positive_product(selected, indices, is_row=True):
    product = 1
    for idx in indices:
        r, c = (selected[idx], idx) if is_row else (idx, selected[idx])
        product *= grid[r][c]
    return product > 0

# Try to find a valid selection
for combination in itertools.combinations(range(36), 11):
    # Convert linear indices to 2D positions
    positions = [(idx // 6, idx % 6) for idx in combination]
    positions.append((0, 5))  # Include the given position

    # Check row and column products
    valid = True
    for i in range(6):
        row_indices = [j for j, (r, c) in enumerate(positions) if r == i]
        col_indices = [j for j, (r, c) in enumerate(positions) if c == i]
        if not is_positive_product(row_indices, row_indices) or not is_positive_product(col_indices, col_indices):
            valid = False
            break

    if valid:
        selected_positions = positions
        break

# Format the output as required
output = ', '.join(f"{r} {c}" for r, c in selected_positions)
print(f"<<<{output}>>>")