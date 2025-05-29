# The grid
grid = [
    [5, -2, 2, -8, -6],
    [1, -2, 7, -1, 8],
    [4, -7, 6, -7, -3],
    [8, -7, 6, 8, 8],
    [-4, -5, 4, 6, 3]
]

# Selected positions from the previous output
selected_positions = [(0, 1), (2, 4), (0, 4), (0, 0), (1, 1), (2, 3), (0, 2), (1, 0), (1, 3)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    # Check rows
    for i in range(5):
        row_product = 1
        for j in range(5):
            if (i, j) in selected_positions:
                row_product *= grid[i][j]
        if row_product <= 0:
            return False
    
    # Check columns
    for j in range(5):
        col_product = 1
        for i in range(5):
            if (i, j) in selected_positions:
                col_product *= grid[i][j]
        if col_product <= 0:
            return False
    
    return True

# Verify the selected positions
is_valid = is_valid_selection(selected_positions)
print(f"Is the selection valid? {is_valid}")