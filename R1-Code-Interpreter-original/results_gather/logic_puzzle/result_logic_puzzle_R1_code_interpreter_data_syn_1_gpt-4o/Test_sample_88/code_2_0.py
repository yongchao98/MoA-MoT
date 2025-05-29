# Define the grid
grid = [
    [-7, 8, 9, 3, 2, -9],
    [-5, -8, 6, 10, 12, 2],
    [-4, 5, 12, 7, 6, 12],
    [1, 7, 7, -2, -5, 6],
    [-10, 6, -5, 9, 2, -8],
    [10, -10, 2, 1, -1, -8]
]

# Initialize the selected positions with the given starting point
selected_positions = [(0, 5)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    # Check if all rows and columns have even negative counts
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Backtracking function to select numbers
def backtrack(r, c):
    if len(selected_positions) == 12:
        return is_valid_selection(selected_positions)
    
    for i in range(r, 6):
        for j in range(c if i == r else 0, 6):
            if (i, j) not in selected_positions:
                selected_positions.append((i, j))
                if backtrack(i, j + 1):
                    return True
                selected_positions.pop()
    return False

# Start backtracking from the initial position
backtrack(0, 0)

# Print the selected positions
print(selected_positions)