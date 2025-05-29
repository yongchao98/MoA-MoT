from itertools import combinations

# Define the grid
grid = [
    [9, 7, 4, 5, 5],
    [-7, 4, -3, -4, 5],
    [7, 1, 3, 2, 8],
    [4, 9, 2, 1, 6],
    [-4, 6, -9, 7, 3]
]

# Pre-select the number at (0,4)
selected_positions = [(0, 4)]

# Function to calculate the product of selected numbers in a row or column
def product_of_selected(selected, indices, is_row=True):
    product = 1
    for idx in indices:
        if is_row:
            product *= grid[selected[idx][0]][selected[idx][1]]
        else:
            product *= grid[selected[idx][1]][selected[idx][0]]
    return product

# Check if the selected positions satisfy the constraints
def is_valid_selection(selected):
    # Check each row
    for row in range(5):
        row_indices = [i for i, (r, c) in enumerate(selected) if r == row]
        if len(row_indices) > 0 and product_of_selected(selected, row_indices) <= 0:
            return False
    
    # Check each column
    for col in range(5):
        col_indices = [i for i, (r, c) in enumerate(selected) if c == col]
        if len(col_indices) > 0 and product_of_selected(selected, col_indices, is_row=False) <= 0:
            return False
    
    return True

# Try to find a valid selection
for comb in combinations([(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)], 9):
    selected = selected_positions + list(comb)
    if is_valid_selection(selected):
        print(selected)
        break