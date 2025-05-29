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
    # Check if white knights are at D3 and D2, black knights at A2 and C1
    return (board[0][3] == 'w' and board[1][3] == 'w' and 
            board[1][0] == 'B' and board[2][2] == 'B')

def find_knights(board, color):
    knights = []
    for i in range(3):  # Fixed range to match board size
        for j in range(4):
            if board[i][j] == color:
                knights.append((i,j))
    return knights

def solve_knight_swap():
    # Initial board (0,0 is top-left, 3x4 board)
    initial_board = [
        ['.', '.', '.', 'B'],  # row 3
        ['w', '.', '.', 'B'],  # row 2
        ['.', '.', 'w', '.']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
        
        current_color = 'w' if is_white_turn else 'B'
        knights = find_knights(board, current_color)
        
        for knight_pos in knights:
            r, c = knight_pos
            valid_moves = get_knight_moves((r,c))
            
            for new_r, new_c in valid_moves:
                if board[new_r][new_c] == '.':  # Empty square
                    new_board = [row[:] for row in board]
                    new_board[new_r][new_c] = current_color
                    new_board[r][c] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Convert to chess notation (A-D for columns, 1-3 for rows)
                        move = f"{current_color},{chr(65+c)}{3-r},{chr(65+new_c)}{3-new_r}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
print("<<<" + ("No" if solution is None else str(solution)) + ">>>")