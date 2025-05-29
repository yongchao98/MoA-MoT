def is_valid_move(start, end, board):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    x, y = ord(position[0]), int(position[1])
    moves = [(x + 2, y + 1), (x + 2, y - 1), (x - 2, y + 1), (x - 2, y - 1),
             (x + 1, y + 2), (x + 1, y - 2), (x - 1, y + 2), (x - 1, y - 2)]
    valid_moves = []
    for mx, my in moves:
        if 'A' <= chr(mx) <= 'D' and 1 <= my <= 3:
            new_pos = chr(mx) + str(my)
            if board.get(new_pos) == '.':
                valid_moves.append(new_pos)
    return valid_moves

def swap_knights(board, moves, turn, depth, max_depth):
    if depth > max_depth:
        return None
    if board['B1'] == 'w' and board['C1'] == 'w' and board['D3'] == 'B' and board['C2'] == 'B':
        return moves

    current_knights = [pos for pos, knight in board.items() if knight == turn]
    for knight_pos in current_knights:
        possible_moves = get_possible_moves(knight_pos, board)
        for move in possible_moves:
            board[move] = turn
            board[knight_pos] = '.'
            moves.append(f"{turn},{knight_pos},{move}")
            next_turn = 'B' if turn == 'w' else 'w'
            result = swap_knights(board, moves, next_turn, depth + 1, max_depth)
            if result:
                return result
            moves.pop()
            board[knight_pos] = turn
            board[move] = '.'
    return None

board = {
    'A3': '.', 'B3': '.', 'C3': '.', 'D3': 'w',
    'A2': '.', 'B2': '.', 'C2': 'w', 'D2': '.',
    'A1': '.', 'B1': 'B', 'C1': 'B', 'D1': '.'
}

max_depth = 20  # Increase the depth limit
moves = swap_knights(board, [], 'w', 0, max_depth)
if moves:
    print(moves)
else:
    print("No")