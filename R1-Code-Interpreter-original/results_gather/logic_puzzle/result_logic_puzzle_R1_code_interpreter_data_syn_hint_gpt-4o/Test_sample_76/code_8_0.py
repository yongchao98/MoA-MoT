# Define the grid
grid = [
    [-8, 10, -8, -3, 2],
    [4, 10, 5, 2, 9],
    [-5, -4, 5, 9, -2],
    [1, 1, 5, -4, -9],
    [-9, -5, 1, 10, 3]
]

# Start with the given number at (0,4)
initial_position = (0, 4)
selected_positions = [initial_position]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 5
    col_neg_count = [0] * 5
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for i in range(5):
        if row_neg_count[i] % 2 != 0 or col_neg_count[i] % 2 != 0:
            return False
    
    return True

# Recursive function to find a valid selection
def find_selection(selected_positions):
    if len(selected_positions) == 10:
        return selected_positions if is_valid_selection(selected_positions) else None
    
    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if is_valid_selection(selected_positions):
                    result = find_selection(selected_positions)
                    if result:
                        return result
                selected_positions.pop()
    return None

# Find the valid selection
result = find_selection(selected_positions)

# Print the selected positions
print(result)