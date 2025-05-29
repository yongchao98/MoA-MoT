# Define the grid
grid = [
    [-8, 5, -9, 6, 7],
    [10, -7, 3, 9, 9],
    [-2, 7, -5, 10, 1],
    [4, 9, 6, -5, 1],
    [7, 9, 3, 7, 5]
]

# Start with the given number 7 at position (0,4)
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

# Try to select 10 numbers
def select_numbers(grid, selected_positions):
    if len(selected_positions) == 10:
        if is_product_positive(selected_positions, grid):
            return selected_positions
        else:
            return None
    
    for i in range(5):
        for j in range(5):
            if (i, j) not in selected_positions:
                selected_positions.append((i, j))
                result = select_numbers(grid, selected_positions)
                if result:
                    return result
                selected_positions.pop()
    return None

# Find the solution
solution = select_numbers(grid, selected_positions)

# Print the solution
print(solution)