import itertools

# The grid
grid = [
    [10, 3, 9, 12, 1, 2],
    [8, 10, -10, -4, 12, -8],
    [6, 12, 12, 9, 11, 10],
    [4, 5, 7, 8, -8, 11],
    [1, -5, 8, -2, 6, 11],
    [6, -11, -3, 11, 9, 3]
]

# Start with the given number at (0,5)
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions):
    # Check rows
    for i in range(6):
        product = 1
        for j in range(6):
            if (i, j) in selected_positions:
                product *= grid[i][j]
        if product <= 0:
            return False

    # Check columns
    for j in range(6):
        product = 1
        for i in range(6):
            if (i, j) in selected_positions:
                product *= grid[i][j]
        if product <= 0:
            return False

    return True

# Try to find a valid selection of 12 numbers
def find_selection():
    # Generate all possible positions except the starting one
    all_positions = [(i, j) for i in range(6) for j in range(6) if (i, j) != (0, 5)]
    
    # Try all combinations of 11 positions (since we already have one)
    for combination in itertools.combinations(all_positions, 11):
        current_selection = selected_positions + list(combination)
        if is_product_positive(current_selection):
            return current_selection

# Find a valid selection
valid_selection = find_selection()

# Format the output
output = ', '.join(f"{i} {j}" for i, j in valid_selection)
print(f"<<<{output}>>>")