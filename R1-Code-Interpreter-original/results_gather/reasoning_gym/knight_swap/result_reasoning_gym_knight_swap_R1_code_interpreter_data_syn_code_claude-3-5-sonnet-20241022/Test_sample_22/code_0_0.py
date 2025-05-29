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

def is_target_reached(board):
    # Check if white knights are at black's initial positions and vice versa
    return (board[2][2] == 'w' and board[2][3] == 'w' and 
            board[0][0] == 'B' and board[1][3] == 'B')

def solve_knight_swap():
    # Initial board state (0-based indexing)
    initial_board = [
        ['w', '.', '.', '.'],
        ['.', '.', 'B', 'w'],
        ['.', '.', '.', 'B']
    ]
    
    # BFS queue: (board, moves_list, is_black_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        # Find all knights of current color
        piece = 'B' if is_black_turn else 'w'
        knights = []
        for r in range(3):
            for c in range(4):
                if board[r][c] == piece:
                    knights.append((r,c))
        
        # Try all possible moves for each knight
        for kr, kc in knights:
            for new_r, new_c in get_knight_moves((kr,kc)):
                if board[new_r][new_c] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[new_r][new_c] = piece
                    new_board[kr][kc] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Convert to chess notation (A1 is bottom-left)
                        from_pos = f"{chr(kc+65)}{3-kr}"
                        to_pos = f"{chr(new_c+65)}{3-new_r}"
                        new_moves = moves + [f"{piece},{from_pos},{to_pos}"]
                        queue.append((new_board, new_moves, not is_black_turn))

    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")