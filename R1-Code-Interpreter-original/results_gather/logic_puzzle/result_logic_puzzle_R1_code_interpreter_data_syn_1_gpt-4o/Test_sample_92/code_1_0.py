# Define the grid
grid = [
    [-4, 2, 3, -10, 11, 2],
    [5, 4, -6, -10, 4, 4],
    [6, -9, 4, 12, 6, 8],
    [12, -10, 6, 6, -5, 1],
    [11, -4, 1, 11, 11, 7],
    [-5, -3, 10, 11, -6, 12]
]

# Start with the given number at (0,5)
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_positive_product(selected_positions):
    # Check rows
    for i in range(6):
        row_product = 1
        for j in range(6):
            if (i, j) in selected_positions:
                row_product *= grid[i][j]
        if row_product <= 0:
            return False
    
    # Check columns
    for j in range(6):
        col_product = 1
        for i in range(6):
            if (i, j) in selected_positions:
                col_product *= grid[i][j]
        if col_product <= 0:
            return False
    
    return True

# Heuristic selection of numbers
def select_numbers():
    # Start with the given number
    selected_positions = [(0, 5)]
    
    # Select positive numbers first
    for i in range(6):
        for j in range(6):
            if (i, j) != (0, 5) and grid[i][j] > 0:
                selected_positions.append((i, j))
                if len(selected_positions) == 12:
                    if is_positive_product(selected_positions):
                        return selected_positions
                    else:
                        selected_positions.pop()
    
    # If not enough positive numbers, balance with negatives
    for i in range(6):
        for j in range(6):
            if (i, j) != (0, 5) and grid[i][j] <= 0:
                selected_positions.append((i, j))
                if len(selected_positions) == 12:
                    if is_positive_product(selected_positions):
                        return selected_positions
                    else:
                        selected_positions.pop()
    
    return selected_positions

# Get the selected positions
selected_positions = select_numbers()

# Print the selected positions
print(selected_positions)