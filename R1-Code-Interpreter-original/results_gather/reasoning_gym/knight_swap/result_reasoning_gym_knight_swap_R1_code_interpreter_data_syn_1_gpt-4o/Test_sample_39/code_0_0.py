def is_valid_move(start, end, board):
    # Calculate the move difference
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if it's a valid L-shaped move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    # Generate all possible L-shaped moves
    x, y = ord(position[0]), int(position[1])
    moves = [(x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
             (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)]
    # Filter valid moves within the board and to empty squares
    valid_moves = []
    for mx, my in moves:
        if 'A' <= chr(mx) <= 'D' and 1 <= my <= 3:
            new_pos = chr(mx) + str(my)
            if board.get(new_pos) == '.':
                valid_moves.append(new_pos)
    return valid_moves

def swap_knights(board, moves, turn):
    # Check if the current board is the target configuration
    target = {'A2': 'w', 'B1': 'w', 'B3': 'w', 'C2': 'B', 'D2': 'B'}
    if all(board[pos] == target[pos] for pos in target):
        return moves

    # Determine whose turn it is
    current_knights = [pos for pos, knight in board.items() if knight.lower() == turn]
    next_turn = 'w' if turn == 'b' else 'b'

    for knight_pos in current_knights:
        possible_moves = get_possible_moves(knight_pos, board)
        for move in possible_moves:
            # Make the move
            board[move] = board[knight_pos]
            board[knight_pos] = '.'
            moves.append(f"{board[move]},{knight_pos},{move}")

            # Recurse
            result = swap_knights(board, moves, next_turn)
            if result:
                return result

            # Undo the move
            moves.pop()
            board[knight_pos] = board[move]
            board[move] = '.'

    return None

# Initial board setup
board = {
    'A2': 'B', 'B1': 'B', 'B3': 'B',
    'C2': 'w', 'D2': 'w',
    'A1': '.', 'A3': '.', 'B2': '.', 'C1': '.', 'C3': '.', 'D1': '.', 'D3': '.'
}

# Attempt to solve the puzzle
moves = swap_knights(board, [], 'b')
if moves:
    print(moves)
else:
    print("No")