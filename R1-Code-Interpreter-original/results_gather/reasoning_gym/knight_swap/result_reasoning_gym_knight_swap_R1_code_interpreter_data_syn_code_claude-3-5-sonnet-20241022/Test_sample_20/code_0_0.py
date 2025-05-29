from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at C2, C3 and black knights at A2, D2
    return (board[1][2] == 'w' and board[0][2] == 'w' and 
            board[1][0] == 'B' and board[1][3] == 'B')

def find_solution():
    # Initial board (0,0 is A1)
    initial_board = [
        ['.', '.', 'B', '.'],  # Row 3
        ['w', '.', 'B', 'w'],  # Row 2
        ['.', '.', '.', '.']   # Row 1
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
        
        visited.add(board_str)
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        # Find all pieces of current color
        pieces = [(r, c) for r in range(3) for c in range(4) 
                 if board[r][c] == piece]
        
        for r, c in pieces:
            for new_r, new_c in get_knight_moves((c, r)):
                if board[new_r][new_c] == '.':  # Empty square
                    new_board = copy.deepcopy(board)
                    new_board[r][c] = '.'
                    new_board[new_r][new_c] = piece
                    move = f"{piece},{chr(65+c)}{3-r},{chr(65+new_c)}{3-new_r}"
                    queue.append((new_board, moves + [move], not is_white_turn))

    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")