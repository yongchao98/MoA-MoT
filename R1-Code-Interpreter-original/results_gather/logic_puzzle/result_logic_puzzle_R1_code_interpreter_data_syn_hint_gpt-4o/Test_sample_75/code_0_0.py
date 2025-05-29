grid = [
    [2, -6, -10, 1, 2],
    [-3, 2, 7, -3, 6],
    [3, 7, 4, 9, 1],
    [-4, 4, 3, 7, 7],
    [8, 4, 3, 5, 10]
]

# Start with the given number
selected_positions = [(0, 4)]
selected_numbers = [2]

# Function to check if adding a number keeps the product positive
def is_valid_selection(selected_positions, new_position):
    row_counts = [0] * 5
    col_counts = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_counts[r] += 1
            col_counts[c] += 1
    r, c = new_position
    if grid[r][c] < 0:
        row_counts[r] += 1
        col_counts[c] += 1
    # Check if any row or column has an odd number of negatives
    for count in row_counts + col_counts:
        if count % 2 != 0:
            return False
    return True

# Try to select 9 more numbers
for r in range(5):
    for c in range(5):
        if (r, c) not in selected_positions:
            if is_valid_selection(selected_positions, (r, c)):
                selected_positions.append((r, c))
                selected_numbers.append(grid[r][c])
                if len(selected_positions) == 10:
                    break
    if len(selected_positions) == 10:
        break

# Output the selected positions
print(selected_positions)