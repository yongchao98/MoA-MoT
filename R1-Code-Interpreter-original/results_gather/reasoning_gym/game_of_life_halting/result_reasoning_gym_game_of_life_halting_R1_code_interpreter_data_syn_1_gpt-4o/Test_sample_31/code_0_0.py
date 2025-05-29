import numpy as np

def get_neighbors_count(board, x, y):
    # Moore neighborhood with wrapping
    neighbors = [
        ((x-1) % board.shape[0], (y-1) % board.shape[1]),
        ((x-1) % board.shape[0], y),
        ((x-1) % board.shape[0], (y+1) % board.shape[1]),
        (x, (y-1) % board.shape[1]),
        (x, (y+1) % board.shape[1]),
        ((x+1) % board.shape[0], (y-1) % board.shape[1]),
        ((x+1) % board.shape[0], y),
        ((x+1) % board.shape[0], (y+1) % board.shape[1])
    ]
    return sum(board[nx, ny] for nx, ny in neighbors)

def game_of_life_halt(initial_board, max_steps):
    board = np.array(initial_board)
    for step in range(max_steps):
        new_board = np.zeros_like(board)
        for x in range(board.shape[0]):
            for y in range(board.shape[1]):
                live_neighbors = get_neighbors_count(board, x, y)
                if board[x, y] == 1:
                    if live_neighbors == 2 or live_neighbors == 3:
                        new_board[x, y] = 1
                else:
                    if live_neighbors == 3:
                        new_board[x, y] = 1
        board = new_board
        if np.sum(board) == 0:
            return True
    return False

initial_board = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

result = game_of_life_halt(initial_board, 20)
print(result)