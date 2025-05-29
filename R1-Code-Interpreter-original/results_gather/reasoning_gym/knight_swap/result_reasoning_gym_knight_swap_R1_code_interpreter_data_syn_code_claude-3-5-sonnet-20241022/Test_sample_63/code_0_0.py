from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(nx,ny) for nx,ny in moves if 0 <= nx < 4 and 0 <= ny < 3]

def pos_to_idx(x, y):
    return y * 4 + x

def get_valid_moves(board, is_white_turn):
    moves = []
    piece = 'w' if is_white_turn else 'B'
    
    for y in range(3):
        for x in range(4):
            if board[pos_to_idx(x,y)] == piece:
                for nx, ny in get_knight_moves((x,y)):
                    if board[pos_to_idx(nx,ny)] == '.':
                        moves.append((x,y,nx,ny))
    return moves

def make_move(board, move):
    x1,y1,x2,y2 = move
    board_list = list(board)
    idx1, idx2 = pos_to_idx(x1,y1), pos_to_idx(x2,y2)
    board_list[idx2] = board_list[idx1]
    board_list[idx1] = '.'
    return ''.join(board_list)

def is_goal(board):
    return board[pos_to_idx(1,2)] == 'B' and board[pos_to_idx(3,2)] == 'B' and \
           board[pos_to_idx(1,1)] == 'w' and board[pos_to_idx(3,1)] == 'w'

def solve_knight_swap():
    initial = '..B.Bw.w...'
    queue = deque([(initial, [], True)])
    seen = {initial}
    
    while queue:
        board, moves, is_white = queue.popleft()
        
        if is_goal(board):
            return moves
            
        for move in get_valid_moves(board, is_white):
            new_board = make_move(board, move)
            if new_board not in seen:
                x1,y1,x2,y2 = move
                move_str = f"{'w' if is_white else 'B'},{chr(65+x1)}{y1+1},{chr(65+x2)}{y2+1}"
                new_moves = moves + [move_str]
                seen.add(new_board)
                queue.append((new_board, new_moves, not is_white))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")