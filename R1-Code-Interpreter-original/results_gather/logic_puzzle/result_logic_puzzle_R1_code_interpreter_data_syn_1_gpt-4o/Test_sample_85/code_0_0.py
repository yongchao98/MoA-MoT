grid = [
    [9, -11, -5, 1, 7, 9],
    [-12, 10, 3, -5, 3, 6],
    [12, -12, -5, 4, 11, 11],
    [4, 12, 8, -6, 7, 7],
    [12, 3, 4, 7, 12, 4],
    [5, 12, 6, 1, 2, 9]
]

selected_positions = [(0, 5)]
selected_numbers = [grid[0][5]]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected, grid, row=None, col=None):
    product = 1
    if row is not None:
        for c in range(6):
            if (row, c) in selected:
                product *= grid[row][c]
    if col is not None:
        for r in range(6):
            if (r, col) in selected:
                product *= grid[r][col]
    return product > 0

# Select numbers ensuring the product of each row and column is positive
for r in range(6):
    for c in range(6):
        if len(selected_positions) >= 12:
            break
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            selected_numbers.append(grid[r][c])
            # Check if the current selection maintains positive product
            if not (is_product_positive(selected_positions, grid, row=r) and is_product_positive(selected_positions, grid, col=c)):
                # If not, remove the last added number
                selected_positions.pop()
                selected_numbers.pop()

# Ensure we have exactly 12 numbers
while len(selected_positions) > 12:
    selected_positions.pop()
    selected_numbers.pop()

print(selected_positions)