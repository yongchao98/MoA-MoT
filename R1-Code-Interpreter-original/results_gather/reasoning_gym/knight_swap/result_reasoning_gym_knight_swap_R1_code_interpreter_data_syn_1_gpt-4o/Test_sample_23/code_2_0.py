from collections import deque

def is_valid_move(start, end, board):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]
    possible_moves = []
    for dx, dy in moves:
        new_x = chr(ord(position[0]) + dx)
        new_y = str(int(position[1]) + dy)
        new_pos = new_x + new_y
        if new_x in 'ABCD' and new_y in '123' and board.get(new_pos) == '.':
            possible_moves.append(new_pos)
    return possible_moves

def swap_knights_bfs(initial_board):
    target_board = {
        'A3': 'B', 'B3': 'B', 'C3': '.', 'D3': '.',
        'A2': '.', 'B2': '.', 'C2': '.', 'D2': 'w',
        'A1': '.', 'B1': 'w', 'C1': '.', 'D1': '.'
    }
    
    queue = deque([(initial_board, [], 'w')])
    visited = set()

    while queue:
        board, moves, turn = queue.popleft()
        board_tuple = tuple(sorted(board.items()))
        
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        if board == target_board:
            return moves

        current_knights = [pos for pos, knight in board.items() if knight == turn]
        next_turn = 'B' if turn == 'w' else 'w'

        for knight_pos in current_knights:
            possible_moves = get_possible_moves(knight_pos, board)
            for move in possible_moves:
                new_board = board.copy()
                new_board[move] = turn
                new_board[knight_pos] = '.'
                new_moves = moves + [f"{turn},{knight_pos},{move}"]
                queue.append((new_board, new_moves, next_turn))

    return None

initial_board = {
    'A3': 'w', 'B3': 'w', 'C3': '.', 'D3': '.',
    'A2': '.', 'B2': '.', 'C2': '.', 'D2': 'B',
    'A1': '.', 'B1': 'B', 'C1': '.', 'D1': '.'
}

moves = swap_knights_bfs(initial_board)

if moves:
    print(moves)
else:
    print("No")