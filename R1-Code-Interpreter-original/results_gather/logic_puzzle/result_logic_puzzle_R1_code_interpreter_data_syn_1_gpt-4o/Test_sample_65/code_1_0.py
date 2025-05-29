# Define the grid
grid = [
    [-4, -5, -7, -6, 6],
    [5, 7, 4, -6, -1],
    [2, -5, -7, 5, 6],
    [4, -6, -3, 3, -3],
    [9, 5, -5, 5, 6]
]

# Initial selection
selected_positions = [(0, 4)]

# Function to check if the product of selected numbers in each row and column is positive
def is_valid_selection(selected_positions):
    row_products = [1] * 5
    col_products = [1] * 5
    for r, c in selected_positions:
        row_products[r] *= grid[r][c]
        col_products[c] *= grid[r][c]
    return all(p > 0 for p in row_products) and all(p > 0 for p in col_products)

# Backtracking function to select numbers
def backtrack(selected_positions):
    if len(selected_positions) == 10:
        if is_valid_selection(selected_positions):
            return selected_positions
        return None
    
    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                result = backtrack(selected_positions)
                if result is not None:
                    return result
                selected_positions.pop()
    return None

# Get the selected positions
selected_positions = backtrack(selected_positions)

# Print the result
print(selected_positions)