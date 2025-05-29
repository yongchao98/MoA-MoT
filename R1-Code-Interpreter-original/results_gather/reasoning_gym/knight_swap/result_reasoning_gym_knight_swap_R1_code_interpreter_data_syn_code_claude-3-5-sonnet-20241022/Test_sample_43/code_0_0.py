from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ','.join(''.join(row) for row in board)

def string_to_board(s):
    rows = s.split(',')
    return [list(row) for row in rows]

def is_target_reached(board):
    return (board[0][3] == 'B' and board[1][3] == 'B' and 
            board[2][1] == 'w' and board[2][2] == 'w')

def find_solution():
    initial_board = [
        ['.', '.', '.', 'w'],
        ['.', '.', '.', 'w'],
        ['.', 'B', 'B', '.']
    ]
    
    # State: (board_string, is_white_turn, moves_list)
    start = (board_to_string(initial_board), True, [])
    visited = set()
    queue = deque([start])
    
    while queue:
        board_str, white_turn, moves = queue.popleft()
        if board_str in visited:
            continue
        visited.add(board_str)
        
        board = string_to_board(board_str)
        if is_target_reached(board):
            return moves
        
        piece = 'w' if white_turn else 'B'
        # Find all pieces of current color
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((j, i))  # (x,y)
        
        # Try all possible moves for each piece
        for px, py in pieces:
            for nx, ny in get_knight_moves((px, py)):
                if board[ny][nx] == '.':  # Empty square
                    new_board = [row[:] for row in board]
                    new_board[py][px] = '.'
                    new_board[ny][nx] = piece
                    
                    move = f"{piece},{chr(65+px)}{3-py},{chr(65+nx)}{3-ny}"
                    new_moves = moves + [move]
                    
                    new_state = (board_to_string(new_board), not white_turn, new_moves)
                    queue.append(new_state)
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")