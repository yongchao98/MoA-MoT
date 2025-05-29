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
    # Adjust valid range for 3x4 board (rows 0-2, cols 0-3)
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # White knights should be at D3, B2 positions (board[0][3] and board[1][1])
    # Black knights should be at A3, A2 positions (board[0][0] and board[1][0])
    return (board[0][3] == 'w' and board[1][1] == 'w' and 
            board[0][0] == 'B' and board[1][0] == 'B')

def get_piece_positions(board, piece):
    return [(i,j) for i in range(3) for j in range(4) if board[i][j] == piece]

def solve_knight_swap():
    # Initial board (3x4)
    initial_board = [
        ['w', '.', '.', 'B'],  # row 3 (index 0)
        ['w', 'B', '.', '.'],  # row 2 (index 1)
        ['.', '.', '.', '.']   # row 1 (index 2)
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
            
        visited.add(board_state)
        
        if is_target_reached(board):
            return moves
            
        current_piece = 'w' if is_white_turn else 'B'
        pieces = get_piece_positions(board, current_piece)
        
        for from_pos in pieces:
            valid_moves = get_knight_moves(from_pos)
            for to_pos in valid_moves:
                if board[to_pos[0]][to_pos[1]] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[to_pos[0]][to_pos[1]] = current_piece
                    new_board[from_pos[0]][from_pos[1]] = '.'
                    
                    # Convert to chess notation (row 0 = row 3, etc.)
                    from_notation = f"{chr(from_pos[1]+65)}{3-from_pos[0]}"
                    to_notation = f"{chr(to_pos[1]+65)}{3-to_pos[0]}"
                    move = f"{current_piece},{from_notation},{to_notation}"
                    
                    queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")