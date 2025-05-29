def move_right(row):
    # Remove zeros and combine numbers
    new_row = [num for num in row if num != 0]
    for i in range(len(new_row) - 1, 0, -1):
        if new_row[i] == new_row[i - 1]:
            new_row[i] *= 2
            new_row[i - 1] = 0
    # Remove zeros again and fill with zeros to the left
    new_row = [num for num in new_row if num != 0]
    return [0] * (len(row) - len(new_row)) + new_row

# Initial grid
grid = [
    [8, 8],
    [8, 32]
]

# First move (r)
grid[0] = move_right(grid[0])
grid[1] = move_right(grid[1])

# Second move (r)
grid[0] = move_right(grid[0])
grid[1] = move_right(grid[1])

print(grid)