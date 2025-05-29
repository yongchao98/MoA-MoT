# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['0', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Place the black piece at (2, 1)
grid[2][0] = '0'

# Check and flip horizontally
for row in range(4):
    for col in range(4):
        if grid[row][col] == '0':
            # Check to the right
            if col + 2 < 4 and grid[row][col + 1] == '1' and grid[row][col + 2] == '0':
                grid[row][col + 1] = '0'
            # Check to the left
            if col - 2 >= 0 and grid[row][col - 1] == '1' and grid[row][col - 2] == '0':
                grid[row][col - 1] = '0'

# Check and flip vertically
for col in range(4):
    for row in range(4):
        if grid[row][col] == '0':
            # Check downwards
            if row + 2 < 4 and grid[row + 1][col] == '1' and grid[row + 2][col] == '0':
                grid[row + 1][col] = '0'
            # Check upwards
            if row - 2 >= 0 and grid[row - 1][col] == '1' and grid[row - 2][col] == '0':
                grid[row - 1][col] = '0'

# Check and flip diagonally
for row in range(4):
    for col in range(4):
        if grid[row][col] == '0':
            # Check diagonal down-right
            if row + 2 < 4 and col + 2 < 4 and grid[row + 1][col + 1] == '1' and grid[row + 2][col + 2] == '0':
                grid[row + 1][col + 1] = '0'
            # Check diagonal down-left
            if row + 2 < 4 and col - 2 >= 0 and grid[row + 1][col - 1] == '1' and grid[row + 2][col - 2] == '0':
                grid[row + 1][col - 1] = '0'
            # Check diagonal up-right
            if row - 2 >= 0 and col + 2 < 4 and grid[row - 1][col + 1] == '1' and grid[row - 2][col + 2] == '0':
                grid[row - 1][col + 1] = '0'
            # Check diagonal up-left
            if row - 2 >= 0 and col - 2 >= 0 and grid[row - 1][col - 1] == '1' and grid[row - 2][col - 2] == '0':
                grid[row - 1][col - 1] = '0'

# Flatten the grid to a single list for the final output
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)