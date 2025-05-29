from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 4]

def board_to_string(board):
    return ','.join(''.join(row) for row in board)

def string_to_board(s):
    rows = s.split(',')
    return [list(row) for row in rows]

def get_positions(board, piece):
    return [(i,j) for i in range(4) for j in range(4) if board[i][j] == piece]

def is_target_reached(board):
    w_pos = set((i,j) for i in range(4) for j in range(4) if board[i][j] == 'w')
    target_w = {(1,0), (3,1)}  # B1, D2 positions
    return w_pos == target_w

def solve():
    # Initial board
    initial = [
        ['.','.','.','.'],
        ['B','.','.','.'],
        ['.','.','.','B'],
        ['w','.','w','.']
    ]
    
    seen = set()
    queue = deque([(board_to_string(initial), [], True)])  # board, moves, white_turn
    
    while queue:
        board_str, moves, white_turn = queue.popleft()
        board = string_to_board(board_str)
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if white_turn else 'B'
        positions = get_positions(board, piece)
        
        for pos in positions:
            for new_pos in get_knight_moves(pos):
                i, j = pos
                ni, nj = new_pos
                if board[ni][nj] == '.':
                    new_board = [row[:] for row in board]
                    new_board[ni][nj] = piece
                    new_board[i][j] = '.'
                    new_board_str = board_to_string(new_board)
                    
                    if new_board_str not in seen:
                        seen.add(new_board_str)
                        move = f"{piece},{chr(ord('A')+j)}{3-i},{chr(ord('A')+nj)}{3-ni}"
                        queue.append((new_board_str, moves + [move], not white_turn))
    
    return None

solution = solve()
if solution:
    print(solution)
else:
    print("No")