from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    x, y = ord(position[0]), int(position[1])
    moves = [(x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
             (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)]
    valid_moves = []
    for mx, my in moves:
        if 'A' <= chr(mx) <= 'D' and 1 <= my <= 3:
            new_pos = chr(mx) + str(my)
            if board.get(new_pos) == '.':
                valid_moves.append(new_pos)
    return valid_moves

def swap_knights_bfs(initial_board):
    target = {'A2': 'w', 'B1': 'w', 'B3': 'w', 'C2': 'B', 'D2': 'B'}
    queue = deque([(initial_board, [], 'b')])
    visited = set()

    while queue:
        board, moves, turn = queue.popleft()
        board_tuple = tuple(sorted(board.items()))
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        if all(board[pos] == target[pos] for pos in target):
            return moves

        current_knights = [pos for pos, knight in board.items() if knight.lower() == turn]
        next_turn = 'w' if turn == 'b' else 'b'

        for knight_pos in current_knights:
            possible_moves = get_possible_moves(knight_pos, board)
            for move in possible_moves:
                new_board = board.copy()
                new_board[move] = new_board[knight_pos]
                new_board[knight_pos] = '.'
                new_moves = moves + [f"{new_board[move]},{knight_pos},{move}"]
                queue.append((new_board, new_moves, next_turn))

    return None

board = {
    'A2': 'B', 'B1': 'B', 'B3': 'B',
    'C2': 'w', 'D2': 'w',
    'A1': '.', 'A3': '.', 'B2': '.', 'C1': '.', 'C3': '.', 'D1': '.', 'D3': '.'
}

moves = swap_knights_bfs(board)
if moves:
    print(moves)
else:
    print("No")