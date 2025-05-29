# Define the grid
grid = [
    [1, 2],
    [3, 1]
]

# Pre-selected number position
selected_positions = [(0, 1)]

# Function to check if the selection is valid
def is_valid_selection(selected_positions):
    row_sums = [0, 0]
    col_sums = [0, 0]
    
    for r, c in selected_positions:
        row_sums[r] += grid[r][c]
        col_sums[c] += grid[r][c]
    
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try selecting two more numbers
for r1 in range(2):
    for c1 in range(2):
        if (r1, c1) not in selected_positions:
            for r2 in range(2):
                for c2 in range(2):
                    if (r2, c2) not in selected_positions and (r2, c2) != (r1, c1):
                        new_selection = selected_positions + [(r1, c1), (r2, c2)]
                        if is_valid_selection(new_selection):
                            print(new_selection)
                            break