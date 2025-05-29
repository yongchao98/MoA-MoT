import numpy as np

def count_neighbors(board, x, y):
    # Get the size of the board
    rows, cols = board.shape
    # Calculate the sum of the 3x3 grid centered on (x, y)
    # Use modulo to wrap around the edges
    return sum(board[(x + i) % rows, (y + j) % cols] for i in range(-1, 2) for j in range(-1, 2)) - board[x, y]

def game_of_life_step(board):
    # Create a copy of the board to store the next state
    new_board = np.copy(board)
    # Iterate over each cell in the board
    for x in range(board.shape[0]):
        for y in range(board.shape[1]):
            # Count the number of live neighbors
            live_neighbors = count_neighbors(board, x, y)
            # Apply the rules of the Game of Life
            if board[x, y] == 1:
                if live_neighbors < 2 or live_neighbors > 3:
                    new_board[x, y] = 0
            else:
                if live_neighbors == 3:
                    new_board[x, y] = 1
    return new_board

# Initial board configuration
initial_board = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,1,0,0,0,0],
    [0,0,0,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,1,0,0,1,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,1,0,0],
    [0,0,0,0,1,0,0,0,0,1],
    [0,0,0,0,0,0,0,0,1,0]
])

# Simulate one step
next_board = game_of_life_step(initial_board)

# Convert the result to a list of lists and print it
print(next_board.tolist())