import itertools

# Define the grid
grid = [
    [10, 6, -4, 3, 2],
    [4, -3, 6, 5, -2],
    [5, -6, 6, -3, -5],
    [7, 5, 6, 9, -6],
    [-4, 6, -5, 4, 8]
]

# Initial selected position
initial_selection = [(0, 4)]

# Function to check if the selection is valid
def is_valid_selection(selection):
    # Check row products
    for i in range(5):
        row_product = 1
        for j in range(5):
            if (i, j) in selection:
                row_product *= grid[i][j]
        if row_product <= 0:
            return False
    
    # Check column products
    for j in range(5):
        col_product = 1
        for i in range(5):
            if (i, j) in selection:
                col_product *= grid[i][j]
        if col_product <= 0:
            return False
    
    return True

# Generate all possible selections of 10 numbers including the initial selection
all_positions = [(i, j) for i in range(5) for j in range(5)]
all_positions.remove((0, 4))  # Remove the initial selection
for selection in itertools.combinations(all_positions, 9):
    full_selection = initial_selection + list(selection)
    if is_valid_selection(full_selection):
        print(full_selection)
        break