def is_valid_move(start, end, board):
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

def swap_knights(board, moves, turn, depth, max_depth):
    target = {'A2': 'w', 'B1': 'w', 'B3': 'w', 'C2': 'B', 'D2': 'B'}
    if all(board[pos] == target[pos] for pos in target):
        return moves

    if depth > max_depth:
        return None

    current_knights = [pos for pos, knight in board.items() if knight.lower() == turn]
    next_turn = 'w' if turn == 'b' else 'b'

    for knight_pos in current_knights:
        possible_moves = get_possible_moves(knight_pos, board)
        for move in possible_moves:
            board[move] = board[knight_pos]
            board[knight_pos] = '.'
            moves.append(f"{board[move]},{knight_pos},{move}")

            result = swap_knights(board, moves, next_turn, depth + 1, max_depth)
            if result:
                return result

            moves.pop()
            board[knight_pos] = board[move]
            board[move] = '.'

    return None

board = {
    'A2': 'B', 'B1': 'B', 'B3': 'B',
    'C2': 'w', 'D2': 'w',
    'A1': '.', 'A3': '.', 'B2': '.', 'C1': '.', 'C3': '.', 'D1': '.', 'D3': '.'
}

max_depth = 20
moves = swap_knights(board, [], 'b', 0, max_depth)
if moves:
    print(moves)
else:
    print("No")