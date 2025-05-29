from itertools import combinations

# Define the grid
grid = [
    [-4, 5, 7, -1],
    [2, -4, -2, -5],
    [6, 6, -4, 4],
    [-1, 2, 1, -1]
]

# Predefined selection
initial_selection = [(0, 3)]

# Function to calculate product of selected numbers in a row or column
def product_of_selected(selected, indices, is_row=True):
    product = 1
    for idx in indices:
        if is_row:
            product *= grid[idx][selected[idx]]
        else:
            product *= grid[selected[idx]][idx]
    return product

# Check if a selection is valid
def is_valid_selection(selection):
    # Check rows
    for i in range(4):
        row_indices = [j for j in range(4) if (i, j) in selection]
        if len(row_indices) > 0 and product_of_selected(row_indices, row_indices) <= 0:
            return False
    
    # Check columns
    for j in range(4):
        col_indices = [i for i in range(4) if (i, j) in selection]
        if len(col_indices) > 0 and product_of_selected(col_indices, col_indices, is_row=False) <= 0:
            return False
    
    return True

# Generate all possible selections of 5 more numbers
all_positions = [(i, j) for i in range(4) for j in range(4) if (i, j) != (0, 3)]
for additional_selection in combinations(all_positions, 5):
    selection = initial_selection + list(additional_selection)
    if is_valid_selection(selection):
        print(selection)
        break