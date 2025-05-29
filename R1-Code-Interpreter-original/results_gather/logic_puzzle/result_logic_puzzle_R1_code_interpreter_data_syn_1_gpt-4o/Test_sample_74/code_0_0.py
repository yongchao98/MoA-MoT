# Define the grid
grid = [
    [-4, -1, -2, 5, 1],
    [4, -7, 6, 5, -3],
    [-6, 10, 9, 3, 10],
    [3, -9, -8, 3, 7],
    [5, 5, 6, 8, -3]
]

# Start with the given number
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    # Check rows
    for i in range(5):
        product = 1
        for j in range(5):
            if (i, j) in selected_positions:
                product *= grid[i][j]
        if product <= 0:
            return False
    
    # Check columns
    for j in range(5):
        product = 1
        for i in range(5):
            if (i, j) in selected_positions:
                product *= grid[i][j]
        if product <= 0:
            return False
    
    return True

# Try to select numbers to satisfy the constraints
def select_numbers(grid, selected_positions):
    # We need to select 10 numbers in total
    while len(selected_positions) < 10:
        for i in range(5):
            for j in range(5):
                if (i, j) not in selected_positions:
                    selected_positions.append((i, j))
                    if is_product_positive(selected_positions, grid):
                        if len(selected_positions) == 10:
                            return selected_positions
                    else:
                        selected_positions.pop()
    return selected_positions

# Get the selected positions
selected_positions = select_numbers(grid, selected_positions)

# Print the selected positions
print(selected_positions)