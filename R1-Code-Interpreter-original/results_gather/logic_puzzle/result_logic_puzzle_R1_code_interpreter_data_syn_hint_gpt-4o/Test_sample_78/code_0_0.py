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

# Try to select 10 numbers
for i in range(5):
    for j in range(5):
        if (i, j) != (0, 4) and len(selected_positions) < 10:
            selected_positions.append((i, j))
            if not is_product_positive(selected_positions, grid):
                selected_positions.pop()

# Print the selected positions
print(selected_positions)