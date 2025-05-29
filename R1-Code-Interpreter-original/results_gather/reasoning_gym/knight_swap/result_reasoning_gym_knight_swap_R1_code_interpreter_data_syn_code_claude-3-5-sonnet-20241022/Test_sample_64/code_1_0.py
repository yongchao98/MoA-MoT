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
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(r,c) for r in range(3) for c in range(4) if board[r][c] == piece]

def is_target_state(board):
    white_pos = get_positions(board, 'w')
    black_pos = get_positions(board, 'B')
    target_white = [(0,0), (0,3)]  # Original black positions
    target_black = [(1,0), (0,1)]  # Original white positions
    return sorted(white_pos) == sorted(target_white) and sorted(black_pos) == sorted(target_black)

def solve_knight_swap():
    # Initial board (3x4)
    initial_board = [
        ['B','w','.','.'],
        ['w','.','.','B'],
        ['.','.','.','.',]
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
            
        current_pieces = get_positions(board, 'w' if is_white_turn else 'B')
        
        for pos in current_pieces:
            for new_pos in get_knight_moves(pos):
                r, c = new_pos
                if board[r][c] == '.':  # Empty square
                    new_board = [row[:] for row in board]
                    old_r, old_c = pos
                    new_board[r][c] = new_board[old_r][old_c]
                    new_board[old_r][old_c] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{'w' if is_white_turn else 'B'},{chr(old_c+65)}{3-old_r},{chr(c+65)}{3-r}"
                        queue.append((new_board, moves + [move], not is_white_turn))

    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")