# Initialize the grid
grid = [
    [0, 16, 0, 0, 0],
    [0, 0, 0, 16, 0],
    [0, 0, 16, 0, 0],
    [0, 0, 32, 0, 0],
    [0, 0, 0, 0, 2]
]

# Initial position of the number 2
position = (4, 4)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'udldu'

# Define the move directions
move_directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Function to perform the move
def move_number(grid, position, move):
    row, col = position
    d_row, d_col = move_directions[move]
    new_row, new_col = row + d_row, col + d_col
    
    # Check if the new position is within bounds
    if 0 <= new_row < 5 and 0 <= new_col < 5:
        # Check if the new position can combine
        if grid[new_row][new_col] == grid[row][col]:
            # Combine the numbers
            grid[new_row][new_col] *= 2
            grid[row][col] = 0
        elif grid[new_row][new_col] == 0:
            # Move the number
            grid[new_row][new_col] = grid[row][col]
            grid[row][col] = 0
        # Update the position
        position = (new_row, new_col)
    
    return grid, position

# Simulate the moves
for move in moves:
    grid, position = move_number(grid, position, move)

# Print the final grid
print(grid)