# Define the grid
grid = [
    [10, 3, 9, 12, 1, 2],
    [8, 10, -10, -4, 12, -8],
    [6, 12, 12, 9, 11, 10],
    [4, 5, 7, 8, -8, 11],
    [1, -5, 8, -2, 6, 11],
    [6, -11, -3, 11, 9, 3]
]

# Start with the given position
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_product_positive(selected_positions, grid):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Recursive function to select numbers
def select_numbers(selected_positions, grid):
    if len(selected_positions) == 12:
        if is_product_positive(selected_positions, grid):
            return selected_positions
        else:
            return None
    
    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if is_product_positive(selected_positions, grid):
                    result = select_numbers(selected_positions, grid)
                    if result is not None:
                        return result
                selected_positions.pop()
    return None

# Find the solution
solution = select_numbers(selected_positions, grid)

# Print the selected positions
print(solution)