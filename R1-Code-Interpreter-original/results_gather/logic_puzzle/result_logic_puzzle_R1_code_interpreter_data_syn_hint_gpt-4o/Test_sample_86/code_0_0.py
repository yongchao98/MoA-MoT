# Define the grid
grid = [
    [-9, 11, 7, 9, 8, 5],
    [1, 3, -3, -3, -4, 12],
    [1, 7, 6, 11, 10, 2],
    [7, 9, 3, 1, 7, 8],
    [11, -3, 5, 12, 1, 9],
    [-3, 8, -1, 6, -5, 5]
]

# Start with the given number at (0,5)
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    # Check rows
    for i in range(6):
        product = 1
        count = 0
        for j in range(6):
            if (i, j) in selected_positions:
                product *= grid[i][j]
                count += 1
        if count > 0 and product <= 0:
            return False
    
    # Check columns
    for j in range(6):
        product = 1
        count = 0
        for i in range(6):
            if (i, j) in selected_positions:
                product *= grid[i][j]
                count += 1
        if count > 0 and product <= 0:
            return False
    
    return True

# Try to select 12 numbers
def select_numbers(grid, selected_positions):
    for i in range(6):
        for j in range(6):
            if (i, j) not in selected_positions:
                selected_positions.append((i, j))
                if len(selected_positions) == 12 and is_product_positive(selected_positions, grid):
                    return selected_positions
                elif len(selected_positions) < 12:
                    result = select_numbers(grid, selected_positions)
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find the solution
solution = select_numbers(grid, selected_positions)

# Print the solution
print(solution)