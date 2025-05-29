# Initial grid setup
grid = [
    [1, 0, '*', '*'],
    ['*', 1, 0, '*'],
    ['*', 0, 1, '*'],
    ['*', '*', '*', '*']
]

# Round 1: Black's move at (1,2)
grid[0][1] = 0
# Flip the piece at (1,1) to black
grid[0][0] = 0

# Round 2: White's move at (1,1)
grid[0][0] = 1
# Flip the piece at (1,2) to white
grid[0][1] = 1

# Convert the grid to a single string
result = ','.join(str(cell) for row in grid for cell in row)
print(result)