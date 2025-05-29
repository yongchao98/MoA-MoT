from collections import deque

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(r,c) for r in range(3) for c in range(4) if board[r][c] == piece]

def is_target_reached(board):
    # Initial white positions: A3, B2
    # Initial black positions: B2, D1
    # Target: whites should be at B2, D1; blacks should be at A3, B2
    w_positions = set((r,c) for r in range(3) for c in range(4) if board[r][c] == 'w')
    b_positions = set((r,c) for r in range(3) for c in range(4) if board[r][c] == 'B')
    
    target_w = {(2,1), (2,3)}  # B2, D1 in array coordinates
    target_b = {(0,0), (1,1)}  # A3, B2 in array coordinates
    
    return w_positions == target_w and b_positions == target_b

def solve_knight_swap():
    # Initialize 3x4 board (0,0 is top-left)
    initial_board = [
        ['w', '.', '.', '.'],  # row 3
        ['B', 'w', '.', '.'],  # row 2
        ['.', '.', '.', 'B']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        current_pieces = get_positions(board, 'w' if is_white_turn else 'B')
        
        for pos in current_pieces:
            for new_pos in get_knight_moves(pos):
                r, c = new_pos
                if board[r][c] == '.':
                    new_board = [row[:] for row in board]
                    new_board[pos[0]][pos[1]] = '.'
                    new_board[r][c] = 'w' if is_white_turn else 'B'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        # Convert to chess notation (A1 is bottom-left)
                        old_pos = f"{chr(pos[1]+65)}{3-pos[0]}"
                        new_pos = f"{chr(c+65)}{3-r}"
                        move = f"{'w' if is_white_turn else 'B'},{old_pos},{new_pos}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return "No"

print(solve_knight_swap())