from collections import deque

def is_valid_move(start, end, board):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    moves = [(2, 1), (2, -1), (-2, 1), (-2, -1),
             (1, 2), (1, -2), (-1, 2), (-1, -2)]
    possible_moves = []
    for dx, dy in moves:
        new_x = chr(ord(position[0]) + dx)
        new_y = str(int(position[1]) + dy)
        new_pos = new_x + new_y
        if 'A' <= new_x <= 'D' and '1' <= new_y <= '3' and board.get(new_pos) == '.':
            possible_moves.append(new_pos)
    return possible_moves

def swap_knights_bfs():
    initial_board = {
        'A1': 'B', 'B1': 'B', 'C1': 'w', 'D1': '.',
        'A2': '.', 'B2': '.', 'C2': '.', 'D2': '.',
        'A3': '.', 'B3': '.', 'C3': 'w', 'D3': '.'
    }
    target = {
        'A1': 'w', 'B1': 'w', 'C1': 'B', 'C3': 'B'
    }
    initial_state = (('C3', 'C1'), ('A1', 'B1'), initial_board, [])
    queue = deque([initial_state])
    visited = set()

    while queue:
        white_knights, black_knights, board, moves = queue.popleft()
        board_tuple = tuple(sorted(board.items()))
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        if all(board[pos] == target[pos] for pos in target):
            return moves

        # White's turn
        for i, pos in enumerate(white_knights):
            possible_moves = get_possible_moves(pos, board)
            for move in possible_moves:
                new_board = board.copy()
                new_board[move] = 'w'
                new_board[pos] = '.'
                new_white_knights = list(white_knights)
                new_white_knights[i] = move
                new_moves = moves + [f"w,{pos},{move}"]
                queue.append((tuple(new_white_knights), black_knights, new_board, new_moves))

        # Black's turn
        for i, pos in enumerate(black_knights):
            possible_moves = get_possible_moves(pos, board)
            for move in possible_moves:
                new_board = board.copy()
                new_board[move] = 'B'
                new_board[pos] = '.'
                new_black_knights = list(black_knights)
                new_black_knights[i] = move
                new_moves = moves + [f"B,{pos},{move}"]
                queue.append((white_knights, tuple(new_black_knights), new_board, new_moves))

    return "No"

result = swap_knights_bfs()
print(result)