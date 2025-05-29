def move_left(grid):
    new_grid = []
    for row in grid:
        # Filter out zeros and move numbers to the left
        new_row = [num for num in row if num != 0]
        # Combine numbers if they are the same
        if len(new_row) > 1 and new_row[0] == new_row[1]:
            new_row = [new_row[0] * 2] + new_row[2:]
        # Fill the rest with zeros
        new_row += [0] * (len(row) - len(new_row))
        new_grid.append(new_row)
    return new_grid

# Initial grid
grid = [
    [0, 4],
    [16, 4]
]

# Perform the moves
grid = move_left(grid)  # First move
grid = move_left(grid)  # Second move

print(grid)