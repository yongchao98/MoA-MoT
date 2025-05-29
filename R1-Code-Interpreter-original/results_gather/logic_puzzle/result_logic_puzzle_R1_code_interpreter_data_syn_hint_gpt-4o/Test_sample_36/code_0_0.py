# Define the grid
grid = [
    [6, -4, 2, 1],
    [-5, 1, -5, -3],
    [-4, 4, -1, 5],
    [-1, 3, 1, -1]
]

# Pre-select the number at (0,3)
selected_positions = [(0, 3)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    # Check rows
    for i in range(4):
        product = 1
        for j in range(4):
            if (i, j) in selected_positions:
                product *= grid[i][j]
        if product <= 0:
            return False
    
    # Check columns
    for j in range(4):
        product = 1
        for i in range(4):
            if (i, j) in selected_positions:
                product *= grid[i][j]
        if product <= 0:
            return False
    
    return True

# Try to select 5 more numbers
def find_selection(grid, selected_positions):
    for i in range(4):
        for j in range(4):
            if (i, j) not in selected_positions:
                selected_positions.append((i, j))
                if len(selected_positions) == 6 and is_product_positive(selected_positions, grid):
                    return selected_positions
                elif len(selected_positions) < 6:
                    result = find_selection(grid, selected_positions)
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find a valid selection
result = find_selection(grid, selected_positions)

# Print the result
if result:
    print("<<<" + ", ".join(f"{i} {j}" for i, j in result) + ">>>")
else:
    print("No valid selection found")