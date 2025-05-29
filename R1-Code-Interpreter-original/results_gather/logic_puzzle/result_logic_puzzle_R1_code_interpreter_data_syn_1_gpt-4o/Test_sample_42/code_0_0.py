import itertools

# Define the grid
grid = [
    [8, 4, 5, 1],
    [1, -4, -1, 4],
    [-5, -3, 7, 6],
    [4, 1, 7, 5]
]

# Initial selected position
selected_positions = [(0, 3)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    # Check rows
    for i in range(4):
        row_product = 1
        for j in range(4):
            if (i, j) in selected_positions:
                row_product *= grid[i][j]
        if row_product <= 0:
            return False
    
    # Check columns
    for j in range(4):
        col_product = 1
        for i in range(4):
            if (i, j) in selected_positions:
                col_product *= grid[i][j]
        if col_product <= 0:
            return False
    
    return True

# Try to find a valid selection
for positions in itertools.combinations([(i, j) for i in range(4) for j in range(4) if (i, j) != (0, 3)], 5):
    current_selection = selected_positions + list(positions)
    if is_valid_selection(current_selection):
        print(current_selection)
        break