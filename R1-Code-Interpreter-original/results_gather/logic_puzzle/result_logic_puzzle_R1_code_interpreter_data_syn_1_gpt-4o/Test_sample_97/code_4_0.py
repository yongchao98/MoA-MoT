# Define the grid
grid = [
    [10, -10, -5, -13, -7, 5, 7],
    [13, 14, 3, 3, 3, -9, 6],
    [2, -10, -9, 5, -3, 12, -13],
    [-3, 14, 12, 9, 11, 6, 5],
    [10, 1, 13, 9, 8, 6, 3],
    [8, 2, -4, 8, 12, 8, 4],
    [13, 1, 14, -2, 1, -13, 1]
]

# Selected positions from the previous output
selected_positions = [(1, 2), (3, 4), (0, 0), (3, 1), (1, 1), (2, 0), (1, 4), (0, 6), (2, 3), (3, 3), (0, 5), (1, 0), (1, 6), (3, 2), (2, 5), (1, 3)]

# Function to check if the product of numbers in a list is positive
def is_product_positive(numbers):
    product = 1
    for num in numbers:
        product *= num
    return product > 0

# Verify the selection
def verify_selection(selected_positions):
    # Check row products
    for row in range(7):
        row_numbers = [grid[row][col] for col in range(7) if (row, col) in selected_positions]
        if len(row_numbers) > 0 and not is_product_positive(row_numbers):
            return False
    
    # Check column products
    for col in range(7):
        col_numbers = [grid[row][col] for row in range(7) if (row, col) in selected_positions]
        if len(col_numbers) > 0 and not is_product_positive(col_numbers):
            return False
    
    return True

# Verify the selected positions
is_valid = verify_selection(selected_positions)

# Print the verification result
print("Selection is valid:", is_valid)