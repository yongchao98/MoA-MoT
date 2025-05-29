def is_valid_move(start, end, board):
    # Calculate the move difference
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    # All possible knight moves
    moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]
    possible_moves = []
    for dx, dy in moves:
        new_x = chr(ord(position[0]) + dx)
        new_y = str(int(position[1]) + dy)
        new_pos = new_x + new_y
        if new_x in 'ABCD' and new_y in '123' and board.get(new_pos) == '.':
            possible_moves.append(new_pos)
    return possible_moves

def swap_knights(board, moves, turn):
    # Check if the current board is the target configuration
    if board['A3'] == 'B' and board['B3'] == 'B' and board['D2'] == 'w' and board['B1'] == 'w':
        return moves

    # Determine whose turn it is
    current_knights = [pos for pos, knight in board.items() if knight == turn]
    next_turn = 'B' if turn == 'w' else 'w'

    for knight_pos in current_knights:
        possible_moves = get_possible_moves(knight_pos, board)
        for move in possible_moves:
            # Make the move
            board[move] = turn
            board[knight_pos] = '.'
            moves.append(f"{turn},{knight_pos},{move}")

            # Recurse
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
    'A3': 'w', 'B3': 'w', 'C3': '.', 'D3': '.',
    'A2': '.', 'B2': '.', 'C2': '.', 'D2': 'B',
    'A1': '.', 'B1': 'B', 'C1': '.', 'D1': '.'
}

# Start the swap process
moves = swap_knights(board, [], 'w')

# Output the result
if moves:
    print(moves)
else:
    print("No")