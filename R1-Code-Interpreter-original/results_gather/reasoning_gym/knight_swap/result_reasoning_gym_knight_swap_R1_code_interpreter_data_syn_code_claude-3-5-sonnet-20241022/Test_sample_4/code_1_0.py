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
    return [(r,c) for r,c in moves if 1 <= r <= 3 and 0 <= c <= 3]  # col from 0 to 3 (A to D)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at B1,B2 and black knights at D1,D2
    return (board[0][1] == '.' and 
            board[2][1] == 'w' and 
            board[1][1] == 'w' and 
            board[2][3] == 'B' and 
            board[1][3] == 'B')

def solve_knight_swap():
    # Initialize board (A=0, B=1, C=2, D=3 for columns)
    initial_board = [
        ['.' for _ in range(4)],  # 0-th row (unused)
        ['.', 'B', '.', 'w'],     # row 1
        ['.', 'B', '.', 'w'],     # row 2
        ['.', '.', '.', '.']      # row 3
    ]
    
    queue = deque([(initial_board, [], 'w')])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, player = queue.popleft()
        
        if is_goal_state(board):
            return moves
            
        # Find all pieces of current player
        pieces = []
        for r in range(1, 4):
            for c in range(4):
                if board[r][c] == player:
                    pieces.append((r, c))
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            r, c = piece_pos
            for new_r, new_c in get_knight_moves(piece_pos):
                if board[new_r][new_c] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[new_r][new_c] = board[r][c]
                    new_board[r][c] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        from_sq = f"{chr(ord('A')+c)}{r}"
                        to_sq = f"{chr(ord('A')+new_c)}{new_r}"
                        new_moves = moves + [f"{player},{from_sq},{to_sq}"]
                        queue.append((new_board, new_moves, 'B' if player == 'w' else 'w'))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("No")
else:
    print(solution)