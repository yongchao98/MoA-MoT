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

def is_target_state(board):
    # Black knights should be where white knights were (A2, B3)
    # White knights should be where black knights were (A3, D2)
    white_positions = set((r,c) for r in range(3) for c in range(4) if board[r][c] == 'w')
    black_positions = set((r,c) for r in range(3) for c in range(4) if board[r][c] == 'B')
    return white_positions == {(0,0), (1,3)} and black_positions == {(1,0), (0,1)}

def solve_knight_swap():
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
        
        for old_r, old_c in current_pieces:
            for new_r, new_c in get_knight_moves((old_r, old_c)):
                if board[new_r][new_c] == '.':
                    new_board = [row[:] for row in board]
                    new_board[new_r][new_c] = new_board[old_r][old_c]
                    new_board[old_r][old_c] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{'w' if is_white_turn else 'B'},{chr(old_c+65)}{3-old_r},{chr(new_c+65)}{3-new_r}"
                        if len(moves) < 15:  # Limit search depth
                            queue.append((new_board, moves + [move], not is_white_turn))

    return "No"

print(solve_knight_swap())