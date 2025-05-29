import numpy as np

def get_neighbors_count(board, x, y):
    # Get the size of the board
    rows, cols = board.shape
    # Calculate the sum of the Moore neighborhood
    # Use modulo for wrapping topology
    neighbors = sum([board[(x+i)%rows, (y+j)%cols] for i in [-1, 0, 1] for j in [-1, 0, 1] if not (i == 0 and j == 0)])
    return neighbors

def game_of_life(board, steps):
    for _ in range(steps):
        new_board = board.copy()
        for x in range(board.shape[0]):
            for y in range(board.shape[1]):
                neighbors = get_neighbors_count(board, x, y)
                if board[x, y] == 1:
                    if neighbors < 2 or neighbors > 3:
                        new_board[x, y] = 0
                else:
                    if neighbors == 3:
                        new_board[x, y] = 1
        board = new_board
        if np.all(board == 0):
            return True
    return False

# Initial board
initial_board = np.array([
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
])

# Check if the game halts within 20 steps
result = game_of_life(initial_board, 20)
print(result)