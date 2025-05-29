import numpy as np

def get_neighbors_count(board, x, y):
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

def game_of_life_halt(initial_board, max_steps=20):
    board = np.array(initial_board)
    for step in range(max_steps):
        new_board = np.zeros_like(board)
        for x in range(board.shape[0]):
            for y in range(board.shape[1]):
                live_neighbors = get_neighbors_count(board, x, y)
                if board[x, y] == 1 and live_neighbors in (2, 3):
                    new_board[x, y] = 1
                elif board[x, y] == 0 and live_neighbors == 3:
                    new_board[x, y] = 1
        if np.all(new_board == 0):
            return True
        board = new_board
    return False

initial_board = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

result = game_of_life_halt(initial_board)
print(result)