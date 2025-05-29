# Initial grid setup
grid = [
    [1, 0, '*', '*'],
    ['*', 1, 0, '*'],
    ['*', 0, 1, '*'],
    ['*', '*', '*', '*']
]

# Round 1: Your move (black) at (1,2)
grid[0][1] = 0
# Flip the piece at (1,0) since it's between two black pieces
grid[0][0] = 0

# Round 2: Opponent's move (white) at (1,1)
grid[0][0] = 1
# Flip the piece at (1,2) since it's between two white pieces
grid[0][1] = 1

# Convert the grid to a single string
result = ','.join(str(cell) for row in grid for cell in row)
print(result)