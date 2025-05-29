grid = [
    [5, 1, 10, 11, 1, 1],
    [7, 9, -1, 12, -9, 8],
    [3, 10, 1, -1, -11, 1],
    [3, -10, -4, -9, 3, 6],
    [3, -2, 3, -10, 2, 3],
    [-9, -3, 8, 4, 6, 6]
]

selected_positions = [(0, 5)]
selected_numbers = [grid[0][5]]

# We need to select 11 more numbers
# Start by selecting positive numbers to ensure positive product
for i in range(6):
    for j in range(6):
        if len(selected_positions) == 12:
            break
        if (i, j) not in selected_positions:
            if grid[i][j] > 0:
                selected_positions.append((i, j))
                selected_numbers.append(grid[i][j])

# If we still need more numbers, ensure even number of negatives in rows/columns
if len(selected_positions) < 12:
    for i in range(6):
        for j in range(6):
            if len(selected_positions) == 12:
                break
            if (i, j) not in selected_positions:
                selected_positions.append((i, j))
                selected_numbers.append(grid[i][j])

# Output the selected positions
print(selected_positions)