# Define the grid
grid = [
    [6, 1, 5, -9, 3, -10],
    [-9, -4, -9, 2, 7, 11],
    [9, -8, 6, 12, 12, -2],
    [9, -2, 11, -6, 8, 12],
    [4, 1, 2, -2, -5, 10],
    [3, -5, 6, -1, 9, 12]
]

# Initial selection
selected_positions = [(0, 5)]
selected_numbers = [-10]

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

# Backtracking function to find a valid selection
def find_selection(selected_positions, grid):
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
                    result = find_selection(selected_positions, grid)
                    if result is not None:
                        return result
                selected_positions.pop()
    return None

# Find a valid selection
result = find_selection(selected_positions, grid)

# Output the selected positions
if result:
    print(', '.join(f"{r} {c}" for r, c in result))
else:
    print("No valid selection found")