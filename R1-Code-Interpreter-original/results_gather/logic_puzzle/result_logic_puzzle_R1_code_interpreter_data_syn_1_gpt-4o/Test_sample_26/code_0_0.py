from itertools import combinations

# Define the grid
grid = [
    [1, 1, 2],
    [3, 4, 3],
    [4, 4, 3]
]

# Pre-selected number position
pre_selected = (0, 2)

# All possible positions except the pre-selected one
positions = [(i, j) for i in range(3) for j in range(3) if (i, j) != pre_selected]

# Function to check if a selection is valid
def is_valid_selection(selected_positions):
    row_sums = [0, 0, 0]
    col_sums = [0, 0, 0]
    
    # Include the pre-selected number
    row_sums[pre_selected[0]] += grid[pre_selected[0]][pre_selected[1]]
    col_sums[pre_selected[1]] += grid[pre_selected[0]][pre_selected[1]]
    
    # Calculate sums for the selected positions
    for pos in selected_positions:
        row_sums[pos[0]] += grid[pos[0]][pos[1]]
        col_sums[pos[1]] += grid[pos[0]][pos[1]]
    
    # Check if all sums are <= 4
    return all(s <= 4 for s in row_sums) and all(s <= 4 for s in col_sums)

# Try all combinations of 3 positions
for selected_positions in combinations(positions, 3):
    if is_valid_selection(selected_positions):
        # Include the pre-selected position in the result
        result = [pre_selected] + list(selected_positions)
        # Format the result as required
        formatted_result = ', '.join(f"{r} {c}" for r, c in result)
        print(f"<<<{formatted_result}>>>")
        break