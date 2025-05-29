# Define the grid
grid = [
    [6, -8, -3, 9, 8],
    [4, 8, 9, 3, -8],
    [9, 4, 5, 2, 1],
    [-4, -2, -4, 2, 7],
    [1, -3, 9, 1, -10]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 4)]

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

# Backtracking function to find a valid selection
def find_selection(selected_positions):
    if len(selected_positions) == 10:
        if is_valid_selection(selected_positions):
            return selected_positions
        else:
            return None
    
    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if is_valid_selection(selected_positions):
                    result = find_selection(selected_positions)
                    if result is not None:
                        return result
                selected_positions.pop()
    
    return None

# Find a valid selection
result = find_selection(selected_positions)

# Print the selected positions
print(result)