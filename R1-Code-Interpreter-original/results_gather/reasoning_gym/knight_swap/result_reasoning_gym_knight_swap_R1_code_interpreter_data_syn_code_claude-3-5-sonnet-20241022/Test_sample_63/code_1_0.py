from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(nx,ny) for nx,ny in moves if 0 <= nx < 4 and 0 <= ny < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def string_to_board(s):
    return [list(s[i:i+4]) for i in range(0, 12, 4)]

def get_valid_moves(board, is_white_turn):
    moves = []
    piece = 'w' if is_white_turn else 'B'
    
    for y in range(3):
        for x in range(4):
            if board[y][x] == piece:
                for nx, ny in get_knight_moves((x,y)):
                    if board[ny][nx] == '.':
                        moves.append((x,y,nx,ny))
    return moves

def make_move(board, move):
    x1,y1,x2,y2 = move
    new_board = [row[:] for row in board]
    new_board[y2][x2] = new_board[y1][x1]
    new_board[y1][x1] = '.'
    return new_board

def is_goal(board):
    return (board[1][1] == 'B' and board[1][3] == 'B' and 
            board[2][1] == 'w' and board[2][3] == 'w')

def solve_knight_swap():
    initial_board = [
        ['.', 'B', '.', 'B'],
        ['.', 'w', '.', 'w'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white = queue.popleft()
        
        if is_goal(board):
            return moves
            
        for move in get_valid_moves(board, is_white):
            new_board = make_move(board, move)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                x1,y1,x2,y2 = move
                move_str = f"{'w' if is_white else 'B'},{chr(65+x1)}{3-y1},{chr(65+x2)}{3-y2}"
                new_moves = moves + [move_str]
                seen.add(board_str)
                queue.append((new_board, new_moves, not is_white))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")