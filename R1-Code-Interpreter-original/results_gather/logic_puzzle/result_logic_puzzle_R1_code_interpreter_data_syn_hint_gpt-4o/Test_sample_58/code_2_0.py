# Define the grid
grid = [
    [-2, 3, 7, 8, 8],
    [-5, 1, 5, 8, 6],
    [7, -3, -1, -4, 2],
    [5, -4, -7, 3, 3],
    [2, -4, 3, 1, 6]
]

# Predefined selected position
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Try to find a valid selection
def find_valid_selection():
    # Start with the given number at (0,4)
    selected_positions = [(0, 4)]
    
    # Manually select numbers to balance the grid
    # Select numbers to ensure positive product in each row and column
    selected_positions.extend([(0, 0), (0, 1), (1, 1), (1, 2), (2, 0), (2, 4), (3, 3), (4, 4)])
    
    if is_valid_selection(selected_positions):
        return selected_positions
    return None

result = find_valid_selection()
print(result)