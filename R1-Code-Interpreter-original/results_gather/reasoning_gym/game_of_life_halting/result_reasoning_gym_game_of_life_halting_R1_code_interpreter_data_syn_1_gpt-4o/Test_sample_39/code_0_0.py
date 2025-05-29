import numpy as np

def count_neighbors(board, x, y):
    # Count live neighbors with wrapping
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

def game_of_life(board, steps):
    for _ in range(steps):
        new_board = np.copy(board)
        for x in range(board.shape[0]):
            for y in range(board.shape[1]):
                live_neighbors = count_neighbors(board, x, y)
                if board[x, y] == 1:
                    if live_neighbors < 2 or live_neighbors > 3:
                        new_board[x, y] = 0
                else:
                    if live_neighbors == 3:
                        new_board[x, y] = 1
        board = new_board
        if np.sum(board) == 0:
            return True
    return False

initial_board = np.array([
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
])

result = game_of_life(initial_board, 20)
print(result)