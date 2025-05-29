# Initialize the grid
grid = [
    [0, 0, 0, 0, 0],
    [16, 0, 8, 32, 0],
    [0, 8, 32, 0, 0],
    [16, 8, 0, 0, 0],
    [0, 32, 0, 2, 0]
]

# Initial position of the number to move
position = (4, 2)

# Sequence of moves
moves = 'llrul'

# Define the movement directions
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Function to perform the move
def move_and_merge(grid, position, move):
    x, y = position
    dx, dy = directions[move]
    new_x, new_y = x + dx, y + dy
    
    # Check if the new position is within bounds
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        # Check if the numbers can merge
        if grid[new_x][new_y] == grid[x][y]:
            grid[new_x][new_y] *= 2
            grid[x][y] = 0
            position = (new_x, new_y)
        elif grid[new_x][new_y] == 0:
            grid[new_x][new_y] = grid[x][y]
            grid[x][y] = 0
            position = (new_x, new_y)
    
    return grid, position

# Simulate the sequence of moves
for move in moves:
    grid, position = move_and_merge(grid, position, move)

# Print the final grid
print(grid)