# Define the grid
grid = [
    [2, 1],
    [1, 3]
]

# Pre-selected number
selected_positions = [(0, 1)]  # (0,1) is already selected

# Check possible selections
def is_valid_selection(selected_positions):
    # Calculate sums for rows and columns
    row_sums = [0, 0]
    col_sums = [0, 0]
    
    for r, c in selected_positions:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    # Check if all sums are <= 4
    return all(s <= 4 for s in row_sums + col_sums)

# Try selecting (0,0) and (1,0)
selected_positions.append((0, 0))
selected_positions.append((1, 0))

# Check if this selection is valid
if is_valid_selection(selected_positions):
    print(selected_positions)
else:
    print("No valid selection found")