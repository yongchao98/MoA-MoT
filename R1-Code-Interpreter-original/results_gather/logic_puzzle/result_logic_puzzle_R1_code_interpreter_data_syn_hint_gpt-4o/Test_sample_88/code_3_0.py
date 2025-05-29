# Define the grid
grid = [
    [-7, 8, 9, 3, 2, -9],
    [-5, -8, 6, 10, 12, 2],
    [-4, 5, 12, 7, 6, 12],
    [1, 7, 7, -2, -5, 6],
    [-10, 6, -5, 9, 2, -8],
    [10, -10, 2, 1, -1, -8]
]

# Initialize the selected positions with the given number
selected_positions = [(0, 5)]

# Function to check if a selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    # Check if all rows and columns have even number of negatives
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Backtracking function to select numbers
def select_numbers(selected_positions, start_row, start_col):
    if len(selected_positions) == 12:
        return selected_positions if is_valid_selection(selected_positions) else None
    
    for r in range(start_row, 6):
        for c in range(start_col if r == start_row else 0, 6):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                result = select_numbers(selected_positions, r, c + 1)
                if result:
                    return result
                selected_positions.pop()
    return None

# Start the selection process
result = select_numbers(selected_positions, 0, 0)

# Print the selected positions
print(result)