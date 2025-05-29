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

def is_target_state(board):
    # Black knights should be where white knights were: (1,3), (2,2)
    # White knights should be where black knights were: (2,3), (1,0)
    return (board[1][3] == 'B' and 
            board[2][2] == 'B' and 
            board[2][3] == 'w' and 
            board[1][0] == 'w')

def find_solution():
    # Initial board state (0-based indexing)
    # Row 3 becomes index 0, Row 1 becomes index 2
    initial_board = [
        ['.','.','.','.'],
        ['B','.','.','.'],
        ['.','.','w','B']
    ]
    
    # Keep track of visited states and moves
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, white_turn
    
    while queue:
        board, moves, white_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if is_target_state(board):
            return moves
            
        # Find all pieces that can move
        piece = 'w' if white_turn else 'B'
        for r in range(3):
            for c in range(4):
                if board[r][c] == piece:
                    # Try all possible moves for this piece
                    for new_r, new_c in get_knight_moves((r,c)):
                        if board[new_r][new_c] == '.':
                            # Make move
                            new_board = copy.deepcopy(board)
                            new_board[r][c] = '.'
                            new_board[new_r][new_c] = piece
                            
                            # Record move
                            # Convert to chess notation (A1 starts from bottom left)
                            from_pos = f"{chr(65+c)}{3-r}"
                            to_pos = f"{chr(65+new_c)}{3-new_r}"
                            new_moves = moves + [f"{piece.lower()},{from_pos},{to_pos}"]
                            
                            queue.append((new_board, new_moves, not white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")