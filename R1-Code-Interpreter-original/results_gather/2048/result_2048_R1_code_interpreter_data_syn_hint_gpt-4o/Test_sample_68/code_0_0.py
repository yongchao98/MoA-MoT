# Initial grid
grid = [
    [0, 16, 0, 0, 0],
    [0, 0, 0, 16, 0],
    [0, 0, 16, 0, 0],
    [0, 0, 32, 0, 0],
    [0, 0, 0, 0, 2]
]

# Initial position of the number 2
position = (4, 4)

# Sequence of moves
moves = 'udldu'

# Function to move the number in the grid
def move_number(grid, position, moves):
    for move in moves:
        x, y = position
        if move == 'u':
            new_position = (x - 1, y)
        elif move == 'd':
            new_position = (x + 1, y)
        elif move == 'l':
            new_position = (x, y - 1)
        elif move == 'r':
            new_position = (x, y + 1)
        
        # Check if the new position is within bounds
        if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
            # Check if we can combine numbers
            if grid[new_position[0]][new_position[1]] == grid[x][y]:
                grid[new_position[0]][new_position[1]] *= 2
                grid[x][y] = 0
            elif grid[new_position[0]][new_position[1]] == 0:
                grid[new_position[0]][new_position[1]] = grid[x][y]
                grid[x][y] = 0
            position = new_position
    return grid

# Apply the moves
result_grid = move_number(grid, position, moves)

# Print the resulting grid
print(result_grid)