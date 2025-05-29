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

def swap_knights(board, moves, turn, depth=0, max_depth=10):
    if board['A3'] == 'B' and board['B3'] == 'B' and board['D2'] == 'w' and board['B1'] == 'w':
        return moves

    if depth > max_depth:
        return None

    current_knights = [pos for pos, knight in board.items() if knight == turn]
    next_turn = 'B' if turn == 'w' else 'w'

    for knight_pos in current_knights:
        possible_moves = get_possible_moves(knight_pos, board)
        for move in possible_moves:
            board[move] = turn
            board[knight_pos] = '.'
            moves.append(f"{turn},{knight_pos},{move}")

            result = swap_knights(board, moves, next_turn, depth + 1, max_depth)
            if result:
                return result

            moves.pop()
            board[knight_pos] = turn
            board[move] = '.'

    return None

board = {
    'A3': 'w', 'B3': 'w', 'C3': '.', 'D3': '.',
    'A2': '.', 'B2': '.', 'C2': '.', 'D2': 'B',
    'A1': '.', 'B1': 'B', 'C1': '.', 'D1': '.'
}

moves = swap_knights(board, [], 'w')

if moves:
    print(moves)
else:
    print("No")