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
    # Initial positions:
    # Black knights: D3, A2
    # White knights: D2, C1
    # After swap:
    # White knights should be at: D3, A2
    # Black knights should be at: D2, C1
    return (board[0][3] == 'w' and  # D3
            board[1][0] == 'w' and  # A2
            board[1][3] == 'B' and  # D2
            board[2][2] == 'B')     # C1

def find_solution():
    # Initial board state (0-based indexing, top to bottom)
    initial_board = [
        ['.','.','.','.'],
        ['B','.','.','.'],
        ['.','.','w','B']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, white_turn
    max_depth = 20  # Reasonable limit for number of moves
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        if len(moves) >= max_depth:
            continue
            
        board_str = board_to_string(board)
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if is_target_state(board):
            return moves
            
        piece = 'w' if white_turn else 'B'
        for r in range(3):
            for c in range(4):
                if board[r][c] == piece:
                    for new_r, new_c in get_knight_moves((r,c)):
                        if board[new_r][new_c] == '.':
                            new_board = copy.deepcopy(board)
                            new_board[r][c] = '.'
                            new_board[new_r][new_c] = piece
                            
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