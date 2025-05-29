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

# Backtracking function to find valid selections
def find_selection(selected_positions, grid, row, col):
    if len(selected_positions) == 12:
        if is_valid_selection(selected_positions, grid):
            return selected_positions
        return None

    for r in range(row, 6):
        for c in range(col if r == row else 0, 6):
            if (r, c) not in selected_positions:
                selected_positions.append((r, c))
                result = find_selection(selected_positions, grid, r, c + 1)
                if result:
                    return result
                selected_positions.pop()
    return None

# Find the valid selection
result = find_selection(selected_positions, grid, 0, 0)

# Print the selected positions
print(result)