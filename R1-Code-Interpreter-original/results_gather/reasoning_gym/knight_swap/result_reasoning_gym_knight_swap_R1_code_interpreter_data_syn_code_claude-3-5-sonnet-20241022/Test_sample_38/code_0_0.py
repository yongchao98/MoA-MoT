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
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(r, c) for r in range(3) for c in range(4) if board[r][c] == piece]

def is_target_state(board):
    white_pos = set((r, c) for r in range(3) for c in range(4) if board[r][c] == 'w')
    return white_pos == {(1, 0), (0, 2)}  # Target positions for white knights

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.', 'w', '.', 'w'],
        ['B', '.', '.', '.'],
        ['.', '.', 'B', '.']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
        visited.add(board_state)
        
        if is_target_state(board):
            return moves
        
        current_pieces = get_positions(board, 'w' if is_white_turn else 'B')
        for start_pos in current_pieces:
            for end_pos in get_knight_moves(start_pos):
                r, c = end_pos
                if board[r][c] == '.':
                    new_board = [row[:] for row in board]
                    # Make move
                    new_board[start_pos[0]][start_pos[1]] = '.'
                    new_board[r][c] = 'w' if is_white_turn else 'B'
                    
                    # Record move
                    from_sq = f"{chr(65+start_pos[1])}{start_pos[0]+1}"
                    to_sq = f"{chr(65+c)}{r+1}"
                    move = f"{'w' if is_white_turn else 'B'},{from_sq},{to_sq}"
                    
                    queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")