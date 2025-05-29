def is_valid_move(start, end, board):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid L-shape move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, board):
    # All possible knight moves
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

def swap_knights():
    # Initial board setup
    board = {
        'A1': 'B', 'B1': 'B', 'C1': 'w', 'D1': '.',
        'A2': '.', 'B2': '.', 'C2': '.', 'D2': '.',
        'A3': '.', 'B3': '.', 'C3': 'w', 'D3': '.'
    }
    # Target positions
    target = {
        'A1': 'w', 'B1': 'w', 'C1': 'B', 'C3': 'B'
    }
    # Moves list
    moves = []
    # Knights positions
    white_knights = ['C3', 'C1']
    black_knights = ['A1', 'B1']
    
    # Simulate moves
    for _ in range(10):  # Limit to 10 moves for simplicity
        # White's turn
        for i, pos in enumerate(white_knights):
            possible_moves = get_possible_moves(pos, board)
            for move in possible_moves:
                if move in target and target[move] == 'w':
                    # Make the move
                    board[move] = 'w'
                    board[pos] = '.'
                    moves.append(f"w,{pos},{move}")
                    white_knights[i] = move
                    break
        # Black's turn
        for i, pos in enumerate(black_knights):
            possible_moves = get_possible_moves(pos, board)
            for move in possible_moves:
                if move in target and target[move] == 'B':
                    # Make the move
                    board[move] = 'B'
                    board[pos] = '.'
                    moves.append(f"B,{pos},{move}")
                    black_knights[i] = move
                    break
        # Check if all knights are in target positions
        if all(board[pos] == target[pos] for pos in target):
            return moves
    return "No"

result = swap_knights()
print(result)