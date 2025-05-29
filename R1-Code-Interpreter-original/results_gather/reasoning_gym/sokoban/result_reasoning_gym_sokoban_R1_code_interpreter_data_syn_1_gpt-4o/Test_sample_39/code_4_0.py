from collections import deque

def is_valid_move(board, x, y):
    return board[x][y] in ['-', 'X']

def move_player(player_pos, direction):
    x, y = player_pos
    if direction == 'U':
        return (x - 1, y)
    elif direction == 'D':
        return (x + 1, y)
    elif direction == 'L':
        return (x, y - 1)
    elif direction == 'R':
        return (x, y + 1)

def move_box(box_pos, direction):
    x, y = box_pos
    if direction == 'U':
        return (x - 1, y)
    elif direction == 'D':
        return (x + 1, y)
    elif direction == 'L':
        return (x, y - 1)
    elif direction == 'R':
        return (x, y + 1)

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@':
                return False
    return True

def sokoban_solver(initial_board, player_start):
    directions = ['U', 'D', 'L', 'R']
    queue = deque([(initial_board, player_start, "")])
    visited = set()

    while queue:
        board, player_pos, path = queue.popleft()
        if is_solved(board):
            return path

        for direction in directions:
            new_player_pos = move_player(player_pos, direction)
            if is_valid_move(board, *new_player_pos):
                new_board = [list(row) for row in board]
                x, y = new_player_pos
                if new_board[x][y] == '@':
                    new_box_pos = move_box(new_player_pos, direction)
                    if is_valid_move(board, *new_box_pos):
                        new_board[x][y] = '-'
                        bx, by = new_box_pos
                        new_board[bx][by] = '@'
                        new_board[player_pos[0]][player_pos[1]] = '-'
                        new_board[x][y] = '*'
                        new_board_tuple = tuple(tuple(row) for row in new_board)
                        if new_board_tuple not in visited:
                            visited.add(new_board_tuple)
                            queue.append((new_board, new_player_pos, path + direction))
                else:
                    new_board[player_pos[0]][player_pos[1]] = '-'
                    new_board[x][y] = '*'
                    new_board_tuple = tuple(tuple(row) for row in new_board)
                    if new_board_tuple not in visited:
                        visited.add(new_board_tuple)
                        queue.append((new_board, new_player_pos, path + direction))

    return "No solution"

initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '*', '@', '-', '-', '-', 'X', '+'],
    ['+', '@', '$', '+', '-', '@', '-', '+'],
    ['+', 'X', '+', 'X', '-', '-', '-', '+'],
    ['+', '+', '+', 'X', '-', 'X', '-', '+'],
    ['+', '+', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', 'X', '@', '-', '+'],
    ['+', '-', '@', '-', '-', '@', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+']
]

player_start = (1, 1)
solution = sokoban_solver(initial_board, player_start)
print(solution)