import itertools

# Define the grid
grid = [
    [6, -4, 2, 1],
    [-5, 1, -5, -3],
    [-4, 4, -1, 5],
    [-1, 3, 1, -1]
]

# Initial selected position
initial_selection = [(0, 3)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selection):
    # Check rows
    for i in range(4):
        row_product = 1
        for j in range(4):
            if (i, j) in selection:
                row_product *= grid[i][j]
        if row_product <= 0:
            return False
    
    # Check columns
    for j in range(4):
        col_product = 1
        for i in range(4):
            if (i, j) in selection:
                col_product *= grid[i][j]
        if col_product <= 0:
            return False
    
    return True

# Generate all possible selections of 5 more numbers
all_positions = [(i, j) for i in range(4) for j in range(4) if (i, j) != (0, 3)]
for selection in itertools.combinations(all_positions, 5):
    full_selection = initial_selection + list(selection)
    if is_valid_selection(full_selection):
        print(full_selection)
        break