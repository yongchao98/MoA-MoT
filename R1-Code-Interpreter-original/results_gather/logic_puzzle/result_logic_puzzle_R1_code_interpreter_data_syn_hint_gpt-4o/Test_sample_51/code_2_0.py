# Define the grid
grid = [
    [3, -7, 5, 6, 6],
    [1, 8, -2, 6, 7],
    [-6, 1, -7, 7, -6],
    [4, -1, -8, -5, -3],
    [-3, 2, -8, -6, 5]
]

# Start with the given number at (0,4)
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection of 9 numbers
def find_valid_selection():
    # Start with the given number at (0,4)
    selected_positions = [(0, 4)]
    
    # Manually select numbers to ensure constraints are met
    # Select numbers from each row ensuring even number of negatives
    selected_positions.extend([(0, 0), (0, 2), (1, 1), (1, 3), (2, 1), (2, 3), (3, 0), (4, 4)])
    
    if is_valid_selection(selected_positions):
        return selected_positions
    else:
        return None

result = find_valid_selection()
print(result)