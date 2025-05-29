# Initial grid
grid = [
    [0, 4, 8, 0],
    [0, 2, 0, 0],
    [0, 32, 2, 0],
    [2, 4, 0, 0]
]

# Initial position of the number to move
position = (2, 2)

# Sequence of moves
moves = 'rrud'

# Function to perform the moves
def move_number(grid, position, moves):
    x, y = position
    number = grid[x][y]
    grid[x][y] = 0  # Clear the original position

    for move in moves:
        if move == 'r':
            while y < 3 and (grid[x][y + 1] == 0 or grid[x][y + 1] == number):
                if grid[x][y + 1] == number:
                    number *= 2
                    grid[x][y + 1] = 0
                y += 1
        elif move == 'l':
            while y > 0 and (grid[x][y - 1] == 0 or grid[x][y - 1] == number):
                if grid[x][y - 1] == number:
                    number *= 2
                    grid[x][y - 1] = 0
                y -= 1
        elif move == 'u':
            while x > 0 and (grid[x - 1][y] == 0 or grid[x - 1][y] == number):
                if grid[x - 1][y] == number:
                    number *= 2
                    grid[x - 1][y] = 0
                x -= 1
        elif move == 'd':
            while x < 3 and (grid[x + 1][y] == 0 or grid[x + 1][y] == number):
                if grid[x + 1][y] == number:
                    number *= 2
                    grid[x + 1][y] = 0
                x += 1

    grid[x][y] = number  # Place the number in the new position
    return grid

# Perform the moves
final_grid = move_number(grid, position, moves)

# Print the final grid
print(final_grid)