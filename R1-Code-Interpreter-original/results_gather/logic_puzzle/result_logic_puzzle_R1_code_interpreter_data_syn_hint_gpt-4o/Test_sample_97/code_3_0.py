# Define the grid
grid = [
    [10, -10, -5, -13, -7, 5, 7],
    [13, 14, 3, 3, 3, -9, 6],
    [2, -10, -9, 5, -3, 12, -13],
    [-3, 14, 12, 9, 11, 6, 5],
    [10, 1, 13, 9, 8, 6, 3],
    [8, 2, -4, 8, 12, 8, 4],
    [13, 1, 14, -2, 1, -13, 1]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 6)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_count = [0] * 7
    col_count = [0] * 7
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    for count in row_count + col_count:
        if count % 2 != 0:
            return False
    return True

# Function to select numbers ensuring the product condition
def select_numbers(grid, selected_positions):
    # Track the number of negative numbers in each row and column
    row_neg_count = [0] * 7
    col_neg_count = [0] * 7
    
    # Mark the initial selection
    r, c = selected_positions[0]
    if grid[r][c] < 0:
        row_neg_count[r] += 1
        col_neg_count[c] += 1
    
    # Select numbers row by row
    for r in range(7):
        for c in range(7):
            if (r, c) not in selected_positions:
                # Check if adding this number keeps the product positive
                if grid[r][c] < 0:
                    row_neg_count[r] += 1
                    col_neg_count[c] += 1
                selected_positions.append((r, c))
                
                # Check if we have selected 16 numbers
                if len(selected_positions) == 16:
                    if is_product_positive(selected_positions, grid):
                        return selected_positions
                    else:
                        selected_positions.pop()
                        if grid[r][c] < 0:
                            row_neg_count[r] -= 1
                            col_neg_count[c] -= 1
                        continue
                
                # If not yet 16, continue selecting
                if len(selected_positions) < 16:
                    result = select_numbers(grid, selected_positions)
                    if result:
                        return result
                
                # Backtrack if needed
                selected_positions.pop()
                if grid[r][c] < 0:
                    row_neg_count[r] -= 1
                    col_neg_count[c] -= 1
    
    return None

# Find the valid selection of numbers
result = select_numbers(grid, selected_positions)

# Print the result
print(result)