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
initial_position = (0, 5)

# Function to check if the product of selected numbers in a row or column is positive
def is_positive_product(positions, is_row=True):
    for i in range(6):
        product = 1
        for r, c in positions:
            if (is_row and r == i) or (not is_row and c == i):
                product *= grid[r][c]
        if product <= 0:
            return False
    return True

# Try to find a valid selection
for combination in itertools.combinations(range(36), 11):
    # Convert linear indices to 2D positions
    positions = [(idx // 6, idx % 6) for idx in combination]
    positions.append(initial_position)  # Include the given position

    # Check row and column products
    if is_positive_product(positions, is_row=True) and is_positive_product(positions, is_row=False):
        selected_positions = positions
        break

# Format the output as required
output = ', '.join(f"{r} {c}" for r, c in selected_positions)
print(f"<<<{output}>>>")