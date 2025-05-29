# Define the grid
grid = [
    [10, 6, -4, 3, 2],
    [4, -3, 6, 5, -2],
    [5, -6, 6, -3, -5],
    [7, 5, 6, 9, -6],
    [-4, 6, -5, 4, 8]
]

# Start with the given number at (0,4)
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

# Backtracking function to find a valid selection
def find_selection(selected_positions, grid):
    if len(selected_positions) == 10:
        return selected_positions if is_product_positive(selected_positions, grid) else None
    
    for i in range(5):
        for j in range(5):
            if (i, j) not in selected_positions:
                selected_positions.append((i, j))
                if is_product_positive(selected_positions, grid):
                    result = find_selection(selected_positions, grid)
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find a valid selection
result = find_selection(selected_positions, grid)

# Print the selected positions
print(result)