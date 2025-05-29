# Initial grid
grid = [
    [0, 0, 0, 16, 0],
    [0, 0, 0, 16, 0],
    [8, 0, 2, 0, 0],
    [0, 0, 0, 0, 0],
    [2, 0, 8, 32, 0]
]

# Initial position of the number to move
position = (4, 2)  # (row, column) in 0-indexed
number = grid[position[0]][position[1]]

# Movement sequence
moves = 'llruu'

# Function to move the number in the grid
def move_number(grid, position, number, moves):
    for move in moves:
        if move == 'l':
            new_position = (position[0], position[1] - 1)
        elif move == 'r':
            new_position = (position[0], position[1] + 1)
        elif move == 'u':
            new_position = (position[0] - 1, position[1])
        elif move == 'd':
            new_position = (position[0] + 1, position[1])
        
        # Check if the new position is within bounds
        if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
            # Check if the number can combine
            if grid[new_position[0]][new_position[1]] == number:
                number *= 2
                grid[new_position[0]][new_position[1]] = 0  # Clear the old position
            position = new_position
    
    # Place the number in the final position
    grid[position[0]][position[1]] = number
    return grid

# Move the number according to the sequence
result_grid = move_number(grid, position, number, moves)

# Print the resulting grid
print(result_grid)