from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    # Knight moves in L-shape: 2 squares in one direction and 1 square perpendicular
    moves = [
        (row-2, col+1), (row-2, col-1),  # Up 2, right/left 1
        (row+2, col+1), (row+2, col-1),  # Down 2, right/left 1
        (row-1, col+2), (row+1, col+2),  # Right 2, up/down 1
        (row-1, col-2), (row+1, col-2)   # Left 2, up/down 1
    ]
    # Filter valid moves within board boundaries
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(str(cell) for cell in row) for row in board)

def is_target_state(board):
    # Check if white knights are at original black positions and vice versa
    white_pos = set()
    black_pos = set()
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                white_pos.add((i,j))
            elif board[i][j] == 'B':
                black_pos.add((i,j))
    
    # Original positions
    orig_white = {(0,2), (2,1)}  # C3, B1
    orig_black = {(0,0), (2,2)}  # B3, C1
    
    return white_pos == orig_black and black_pos == orig_white

def pos_to_algebraic(row, col):
    # Convert array indices to algebraic notation
    # row 0,1,2 -> 3,2,1
    # col 0,1,2,3 -> A,B,C,D
    return f"{chr(ord('A') + col)}{3-row}"

def solve_knights():
    # Initial board state
    initial = [
        ['B', '.', 'w', ' '],
        ['.', ' ', '.', ' '],
        [' ', 'w', 'B', '.']
    ]
    
    # BFS queue: (board_state, moves_list, white_turn)
    queue = deque([(initial, [], True)])
    seen = {board_to_string(initial)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        # Check if we reached target state
        if is_target_state(board):
            return moves
            
        # Current player's pieces
        piece = 'w' if white_turn else 'B'
        
        # Find all pieces of current color
        pieces = [(r,c) for r in range(3) for c in range(4) if board[r][c] == piece]
        
        # Try all possible moves for each piece
        for pos in pieces:
            row, col = pos
            for new_row, new_col in get_knight_moves(pos):
                # Check if target square is empty
                if board[new_row][new_col] == '.':
                    # Make the move
                    new_board = copy.deepcopy(board)
                    new_board[row][col] = '.'
                    new_board[new_row][new_col] = piece
                    
                    # Generate move notation
                    move = f"{piece},{pos_to_algebraic(row,col)},{pos_to_algebraic(new_row,new_col)}"
                    
                    # Check if new state was seen before
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        queue.append((new_board, moves + [move], not white_turn))
    
    return None

# Solve the puzzle
solution = solve_knights()
print(solution if solution else "No")