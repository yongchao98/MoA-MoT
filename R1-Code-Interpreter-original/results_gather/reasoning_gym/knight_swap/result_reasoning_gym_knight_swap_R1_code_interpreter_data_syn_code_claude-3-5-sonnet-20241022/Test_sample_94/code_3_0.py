from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row-2, col+1), (row-2, col-1),
        (row+2, col+1), (row+2, col-1),
        (row-1, col+2), (row+1, col+2),
        (row-1, col-2), (row+1, col-2)
    ]
    # Only return moves within board boundaries
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(cell for cell in row) for row in board)

def is_target_state(board):
    # Original positions of white knights: C3, B1
    # Original positions of black knights: B3, C1
    # Need to swap these positions
    white_pos = set()
    black_pos = set()
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                white_pos.add((i,j))
            elif board[i][j] == 'B':
                black_pos.add((i,j))
    
    # Target positions (swapped)
    target_white = {(0,1), (2,2)}  # B3, C1
    target_black = {(0,2), (2,1)}  # C3, B1
    
    return white_pos == target_white and black_pos == target_black

def pos_to_algebraic(row, col):
    return f"{chr(ord('A') + col)}{3-row}"

def solve_knights():
    # Initial board state (3x4)
    initial = [
        ['B', '.', 'w', ' '],  # Row 3: B3, empty, C3, blank
        ['.', ' ', '.', ' '],  # Row 2: empty, blank, empty, blank
        [' ', 'w', 'B', '.']   # Row 1: blank, B1, C1, empty
    ]
    
    queue = deque([(initial, [], True)])  # (board, moves, white_turn)
    seen = {board_to_string(initial)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
        
        piece = 'w' if white_turn else 'B'
        
        # Find all current player's pieces
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for row, col in pieces:
            for new_row, new_col in get_knight_moves((row, col)):
                # Check if target square is empty (either '.' or ' ')
                if board[new_row][new_col] in ['.', ' ']:
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[row][col] = '.' if (row, col) in [(0,1), (0,2), (2,1), (2,2)] else ' '
                    new_board[new_row][new_col] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece},{pos_to_algebraic(row,col)},{pos_to_algebraic(new_row,new_col)}"
                        queue.append((new_board, moves + [move], not white_turn))
    
    return None

# Solve the puzzle
solution = solve_knights()
print(solution if solution else "No")