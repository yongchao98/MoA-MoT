import numpy as np

def count_neighbors(board, x, y):
    # Count live neighbors in a wrapping topology
    neighbors = [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1),         (0, 1),
        (1, -1), (1, 0), (1, 1)
    ]
    count = 0
    for dx, dy in neighbors:
        nx, ny = (x + dx) % board.shape[0], (y + dy) % board.shape[1]
        count += board[nx, ny]
    return count

def game_of_life_step(board):
    new_board = np.zeros_like(board)
    for x in range(board.shape[0]):
        for y in range(board.shape[1]):
            live_neighbors = count_neighbors(board, x, y)
            if board[x, y] == 1:
                if live_neighbors == 2 or live_neighbors == 3:
                    new_board[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_board[x, y] = 1
    return new_board

def will_halt(initial_board, max_steps=20):
    board = np.array(initial_board)
    for step in range(max_steps):
        board = game_of_life_step(board)
        if np.all(board == 0):
            return True
    return False

initial_board = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

result = will_halt(initial_board)
print(result)