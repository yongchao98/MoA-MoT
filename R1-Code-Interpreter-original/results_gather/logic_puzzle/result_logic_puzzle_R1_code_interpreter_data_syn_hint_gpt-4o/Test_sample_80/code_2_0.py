# Define the grid
grid = [
    [1, -5, -7, 1, -2, 1],
    [-9, 3, -4, -3, 9, 6],
    [2, 12, -1, -9, 5, 11],
    [10, 11, 11, 10, 2, 12],
    [2, 9, 4, 6, 3, 7],
    [1, -5, 5, 6, 1, 1]
]

# Start with the given position
selected_positions = [(0, 5)]

# Function to check if the product of selected numbers in a row or column is positive
def is_valid_selection(selected_positions, grid):
    row_count = [0] * 6
    col_count = [0] * 6
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    for count in row_count + col_count:
        if count % 2 != 0:
            return False
    return True

# Greedy selection of numbers
def select_numbers(grid):
    selected_positions = [(0, 5)]
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6

    # Mark the initial position
    if grid[0][5] < 0:
        row_neg_count[0] += 1
        col_neg_count[5] += 1

    for r in range(6):
        for c in range(6):
            if (r, c) not in selected_positions:
                if len(selected_positions) < 12:
                    selected_positions.append((r, c))
                    if grid[r][c] < 0:
                        row_neg_count[r] += 1
                        col_neg_count[c] += 1

    # Ensure the selection is valid
    if is_valid_selection(selected_positions, grid):
        return selected_positions
    else:
        return None

# Find the valid selection
result = select_numbers(grid)

# Print the selected positions
print(result)