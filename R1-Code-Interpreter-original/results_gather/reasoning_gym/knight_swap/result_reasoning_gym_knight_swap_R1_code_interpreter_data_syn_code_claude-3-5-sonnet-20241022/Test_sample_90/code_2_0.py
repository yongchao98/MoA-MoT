from collections import deque
import copy

def get_knight_moves(pos):
    moves = [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]
    row, col = pos
    valid = []
    for dx, dy in moves:
        new_row = row + dx
        new_col = col + dy
        if 0 <= new_row < 3 and 0 <= new_col < 4:  # 3x4 board
            valid.append((new_row, new_col))
    return valid

def pos_to_notation(pos):
    row, col = pos
    return chr(col + ord('A')) + str(3 - row)  # Convert to chess notation

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def find_solution():
    # Initial board (0,0 is top-left)
    initial = [
        ['B', '.', '.', '.'],  # row 3
        ['B', '.', '.', '.'],  # row 2
        ['.', '.', 'w', 'w']   # row 1
    ]
    
    # Queue for BFS: (board, moves_list, is_white_turn)
    queue = deque([(initial, [], True)])
    visited = {board_to_string(initial)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        # Find current positions of all knights
        knights = []
        for row in range(3):
            for col in range(4):
                if board[row][col] in ['w', 'B']:
                    knights.append((row, col, board[row][col]))
        
        # Try all possible moves for current player
        for row, col, piece in knights:
            if (piece == 'w') == white_turn:  # Current player's pieces
                for new_row, new_col in get_knight_moves((row, col)):
                    if board[new_row][new_col] == '.':  # Empty square
                        # Make move
                        new_board = copy.deepcopy(board)
                        new_board[row][col] = '.'
                        new_board[new_row][new_col] = piece
                        
                        board_str = board_to_string(new_board)
                        if board_str not in visited:
                            visited.add(board_str)
                            
                            # Record move
                            move = f"{piece},{pos_to_notation((row,col))},{pos_to_notation((new_row,new_col))}"
                            new_moves = moves + [move]
                            
                            # Check if target reached (all pieces in bottom row)
                            if new_board[2][0] == 'B' and new_board[2][1] == 'B' and \
                               new_board[2][2] == 'w' and new_board[2][3] == 'w':
                                return new_moves
                            
                            queue.append((new_board, new_moves, not white_turn))
    
    return "No"

solution = find_solution()
print(solution)