from collections import deque

def get_valid_knight_moves(x, y, board):
    # Knight move patterns
    patterns = [
        (-2,-1), (-2,1), (2,-1), (2,1),
        (-1,-2), (-1,2), (1,-2), (1,2)
    ]
    valid_moves = []
    for dx, dy in patterns:
        new_x, new_y = x + dx, y + dy
        # Check bounds (0-3 for x, 0-2 for y)
        if 0 <= new_x < 4 and 0 <= new_y < 3:
            # Check if target square is empty
            if board[new_y][new_x] == '.':
                valid_moves.append((new_x, new_y))
    return valid_moves

def get_piece_positions(board, piece):
    positions = []
    for y in range(3):
        for x in range(4):
            if board[y][x] == piece:
                positions.append((x, y))
    return positions

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def make_move(board, from_pos, to_pos):
    new_board = [row[:] for row in board]
    piece = new_board[from_pos[1]][from_pos[0]]
    new_board[from_pos[1]][from_pos[0]] = '.'
    new_board[to_pos[1]][to_pos[0]] = piece
    return new_board

def is_final_state(board):
    # White knights should be at B3 and C1
    # Black knights should be at B1 and D2
    white_at = set((x, y) for y in range(3) for x in range(4) if board[y][x] == 'w')
    black_at = set((x, y) for y in range(3) for x in range(4) if board[y][x] == 'B')
    return white_at == {(1,0), (2,2)} and black_at == {(1,2), (3,1)}

def solve():
    # Initial board (using correct coordinates)
    initial_board = [
        ['.', 'B', '.', '.'],  # row 3 (index 0)
        ['.', '.', '.', 'w'],  # row 2 (index 1)
        ['.', 'w', 'B', '.']   # row 1 (index 2)
    ]
    
    queue = deque([(initial_board, [], True)])  # (board, moves, is_white_turn)
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, white_turn = queue.popleft()
        
        if is_final_state(current_board):
            return moves
        
        # Get current player's pieces
        current_piece = 'w' if white_turn else 'B'
        piece_positions = get_piece_positions(current_board, current_piece)
        
        # Try all possible moves for current player's pieces
        for pos in piece_positions:
            valid_moves = get_valid_knight_moves(pos[0], pos[1], current_board)
            for new_pos in valid_moves:
                # Make the move
                new_board = make_move(current_board, pos, new_pos)
                board_str = board_to_string(new_board)
                
                if board_str not in visited:
                    visited.add(board_str)
                    # Format move as per requirements
                    move = f"{current_piece},{chr(65+pos[0])}{3-pos[1]},{chr(65+new_pos[0])}{3-new_pos[1]}"
                    queue.append((new_board, moves + [move], not white_turn))
    
    return "No"

result = solve()
print(result)