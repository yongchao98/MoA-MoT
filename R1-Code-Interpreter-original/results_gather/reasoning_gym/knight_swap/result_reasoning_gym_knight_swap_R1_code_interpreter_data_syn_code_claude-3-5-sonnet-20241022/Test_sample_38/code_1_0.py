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
    return [(r, c) for r, c in moves if 0 <= r < 4 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(r, c) for r in range(4) for c in range(4) if board[r][c] == piece]

def is_target_state(board):
    # White knights should be at B2 and C1 positions
    white_positions = set((r, c) for r in range(4) for c in range(4) if board[r][c] == 'w')
    target_positions = {(1, 1), (0, 2)}  # B2 and C1 in zero-based indices
    return white_positions == target_positions

def solve_knight_swap():
    # Initial board (0-based indices)
    initial_board = [
        ['.', 'w', '.', 'w'],  # row 3
        ['B', '.', '.', '.'],  # row 2
        ['.', '.', 'B', '.']   # row 1
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
        visited.add(board_state)
        
        if len(moves) > 20:  # Limit search depth
            continue
            
        if is_target_state(board):
            return moves
        
        current_pieces = get_positions(board, 'w' if is_white_turn else 'B')
        for start_pos in current_pieces:
            for end_pos in get_knight_moves(start_pos):
                r, c = end_pos
                if r < len(board) and c < len(board[0]) and board[r][c] == '.':
                    new_board = [row[:] for row in board]
                    # Make move
                    new_board[start_pos[0]][start_pos[1]] = '.'
                    new_board[r][c] = 'w' if is_white_turn else 'B'
                    
                    # Convert to algebraic notation (A1 is bottom-left)
                    from_sq = f"{chr(65+start_pos[1])}{3-start_pos[0]}"
                    to_sq = f"{chr(65+c)}{3-r}"
                    move = f"{'w' if is_white_turn else 'B'},{from_sq},{to_sq}"
                    
                    queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")