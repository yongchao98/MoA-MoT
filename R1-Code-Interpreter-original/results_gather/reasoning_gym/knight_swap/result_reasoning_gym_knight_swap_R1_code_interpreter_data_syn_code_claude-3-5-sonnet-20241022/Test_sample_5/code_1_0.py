from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    # Knight moves in board coordinates (row first, then column)
    moves = [(row-2,col+1), (row-2,col-1), (row+2,col+1), (row+2,col-1),
             (row-1,col+2), (row-1,col-2), (row+1,col+2), (row+1,col-2)]
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]  # Fixed board dimensions

def pos_to_notation(pos):
    row, col = pos
    return chr(65 + col) + str(3 - row)  # Convert to chess notation

def notation_to_pos(notation):
    col = ord(notation[0]) - 65
    row = 3 - int(notation[1])
    return (row, col)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if black knights are at white starting positions and vice versa
    return (board[2][0] == 'B' and board[2][2] == 'B' and 
            board[1][1] == 'w' and board[1][3] == 'w')

def solve_knight_swap():
    # Initial board state (row-major order)
    initial_board = [
        ['.', '.', '.', '.'],  # row 3
        ['.', 'B', '.', 'B'],  # row 2
        ['w', '.', 'w', '.']   # row 1
    ]
    
    seen = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_black_turn
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in seen:
            continue
        seen.add(board_str)
        
        if is_target_reached(board):
            return moves
            
        piece = 'B' if is_black_turn else 'w'
        
        # Find all pieces of current color
        pieces = [(r,c) for r in range(3) for c in range(4) 
                 if board[r][c] == piece]
        
        # Try all possible moves for each piece
        for pos in pieces:
            for new_pos in get_knight_moves(pos):
                if board[new_pos[0]][new_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[pos[0]][pos[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece
                    
                    # Record move
                    move = f"{piece},{pos_to_notation(pos)},{pos_to_notation(new_pos)}"
                    new_moves = moves + [move]
                    
                    queue.append((new_board, new_moves, not is_black_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")