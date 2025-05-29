# Initial grid
grid = [
    [0, 16, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [32, 2, 8, 32, 0],
    [0, 32, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Function to move and combine numbers in the grid
def move_and_combine(grid, sequence):
    # Define the movement directions
    directions = {
        'r': (0, 1),
        'd': (1, 0),
        'l': (0, -1),
        'u': (-1, 0)
    }
    
    # Function to apply a single move
    def apply_move(grid, direction):
        dx, dy = directions[direction]
        new_grid = [[0] * 5 for _ in range(5)]
        if direction in ['r', 'l']:
            for x in range(5):
                line = [grid[x][y] for y in range(5) if grid[x][y] != 0]
                if direction == 'r':
                    line = line[::-1]
                new_line = []
                skip = False
                for i in range(len(line)):
                    if skip:
                        skip = False
                        continue
                    if i + 1 < len(line) and line[i] == line[i + 1]:
                        new_line.append(line[i] * 2)
                        skip = True
                    else:
                        new_line.append(line[i])
                if direction == 'r':
                    new_line = new_line[::-1]
                for y in range(len(new_line)):
                    new_grid[x][y if direction == 'l' else 4 - y] = new_line[y]
        else:
            for y in range(5):
                line = [grid[x][y] for x in range(5) if grid[x][y] != 0]
                if direction == 'd':
                    line = line[::-1]
                new_line = []
                skip = False
                for i in range(len(line)):
                    if skip:
                        skip = False
                        continue
                    if i + 1 < len(line) and line[i] == line[i + 1]:
                        new_line.append(line[i] * 2)
                        skip = True
                    else:
                        new_line.append(line[i])
                if direction == 'd':
                    new_line = new_line[::-1]
                for x in range(len(new_line)):
                    new_grid[x if direction == 'u' else 4 - x][y] = new_line[x]
        return new_grid
    
    # Apply the sequence of moves
    for move in sequence:
        grid = apply_move(grid, move)
    
    return grid

# Apply the movement and combination
final_grid = move_and_combine(grid, 'rdrrl')

# Print the final grid
print(final_grid)