def is_valid_move(start, end, board):
    # Calculate the move difference
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    # Generate all possible knight moves from a given position
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

def swap_knights(board, moves, turn):
    # Check if the board is in the target configuration
    if board['B1'] == 'w' and board['C1'] == 'w' and board['D3'] == 'B' and board['C2'] == 'B':
        return moves

    # Determine whose turn it is
    current_knights = [pos for pos, knight in board.items() if knight == turn]
    for knight_pos in current_knights:
        possible_moves = get_possible_moves(knight_pos, board)
        for move in possible_moves:
            # Make the move
            board[move] = turn
            board[knight_pos] = '.'
            moves.append(f"{turn},{knight_pos},{move}")
            # Alternate turn
            next_turn = 'B' if turn == 'w' else 'w'
            result = swap_knights(board, moves, next_turn)
            if result:
                return result
            # Undo the move
            moves.pop()
            board[knight_pos] = turn
            board[move] = '.'
    return None

# Initial board setup
board = {
    'A3': '.', 'B3': '.', 'C3': '.', 'D3': 'w',
    'A2': '.', 'B2': '.', 'C2': 'w', 'D2': '.',
    'A1': '.', 'B1': 'B', 'C1': 'B', 'D1': '.'
}

# Start the swap process
moves = swap_knights(board, [], 'w')
if moves:
    print(moves)
else:
    print("No")