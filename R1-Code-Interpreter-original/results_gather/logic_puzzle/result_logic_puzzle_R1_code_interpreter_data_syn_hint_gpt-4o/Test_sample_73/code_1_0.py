# Define the grid
grid = [
    [-8, 2, -9, 9, -6],
    [4, 2, -2, 7, 7],
    [-7, 6, 8, -4, -5],
    [10, -1, -7, 4, 1],
    [-8, -9, 3, 7, 10]
]

# Start with the given number
selected_positions = [(0, 4)]
selected_numbers = [-6]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_neg_count = [0] * 5
    col_neg_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    # Check if all rows and columns have even number of negatives
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

# Recursive function to find a valid selection
def find_selection(selected_positions):
    if len(selected_positions) == 10:
        return is_valid_selection(selected_positions)
    
    for r in range(5):
        for c in range(5):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                if find_selection(selected_positions):
                    return True
                selected_positions.pop()
    return False

# Find a valid selection
find_selection(selected_positions)

# Print the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))