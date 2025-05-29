import itertools

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

# Initial selected position
selected_positions = [(0, 6)]

# Function to check if the product of numbers in a list is positive
def is_product_positive(numbers):
    product = 1
    for num in numbers:
        product *= num
    return product > 0

# Function to check if the current selection satisfies the conditions
def check_selection(selected_positions):
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

# Try to select 16 numbers
def select_numbers():
    # Start with the initial position
    positions = [(r, c) for r in range(7) for c in range(7)]
    positions.remove((0, 6))
    
    # Try all combinations of 15 additional positions
    for additional_positions in itertools.combinations(positions, 15):
        current_selection = selected_positions + list(additional_positions)
        if check_selection(current_selection):
            return current_selection

# Get the selected positions
selected_positions = select_numbers()

# Print the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))