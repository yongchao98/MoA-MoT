import itertools

# Initialize the grid
grid = [
    [13, -5, -8, 5, -11, 11, 8],
    [9, 4, -14, 3, 3, 12, 5],
    [3, 3, -5, 9, 11, 5, 10],
    [6, 10, 3, -2, 7, 14, 6],
    [6, 7, -2, -14, -1, 10, 7],
    [1, -10, 6, 13, 6, -10, 10],
    [-8, -13, 9, 6, 8, 3, 1]
]

# Start with the given number
selected_positions = [(0, 6)]

# Function to check if the product of numbers in a list is positive
def is_product_positive(numbers):
    product = 1
    for number in numbers:
        product *= number
    return product > 0

# Try to select 15 more numbers
for positions in itertools.combinations([(i, j) for i in range(7) for j in range(7) if (i, j) != (0, 6)], 15):
    selected_positions = [(0, 6)] + list(positions)
    
    # Check row and column products
    row_products = [1] * 7
    col_products = [1] * 7
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    if all(p > 0 for p in row_products) and all(p > 0 for p in col_products):
        break

# Output the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))