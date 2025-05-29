from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 0 <= r < 4 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(r,c) for r in range(4) for c in range(4) if board[r][c] == piece]

def is_target_reached(board):
    w_pos = set((r,c) for r in range(4) for c in range(4) if board[r][c] == 'w')
    b_pos = set((r,c) for r in range(4) for c in range(4) if board[r][c] == 'B')
    
    # Target positions (row, col)
    initial_b = {(1,1), (3,0)}  # Original black positions
    initial_w = {(0,0), (1,1)}  # Original white positions
    
    return w_pos == initial_b and b_pos == initial_w

def solve_knight_swap():
    # Initialize board (0,0 is top-left)
    initial_board = [
        ['w', '.', '.', '.'],
        ['.', 'B', '.', '.'],
        ['.', '.', '.', 'B'],
        ['.', '.', '.', '.']
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
                        old_pos = f"{chr(pos[1]+65)}{4-pos[0]}"
                        new_pos = f"{chr(c+65)}{4-r}"
                        move = f"{'w' if is_white_turn else 'B'},{old_pos},{new_pos}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return "No"

print(solve_knight_swap())