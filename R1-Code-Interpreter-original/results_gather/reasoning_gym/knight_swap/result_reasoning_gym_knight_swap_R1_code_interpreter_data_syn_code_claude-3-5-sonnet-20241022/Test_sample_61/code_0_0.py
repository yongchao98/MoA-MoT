from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_piece_positions(board, piece):
    return [(x,y) for y in range(3) for x in range(4) if board[y][x] == piece]

def is_target_state(board):
    white_pos = set(get_piece_positions(board, 'w'))
    return white_pos == {(1,0), (2,0)}  # B1, C1 positions

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.', 'w', '.', '.'],
        ['.', '.', '.', '.'],
        ['.', 'B', 'B', 'w']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], 'w')])  # board, moves, next_color
    
    while queue:
        board, moves, color = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        if is_target_state(board) and len(get_piece_positions(board, 'B')) == 2:
            return moves
            
        pieces = get_piece_positions(board, color)
        for px, py in pieces:
            for nx, ny in get_knight_moves((px, py)):
                if board[ny][nx] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[py][px] = '.'
                    new_board[ny][nx] = color
                    
                    from_pos = chr(px + ord('A')) + str(3-py)
                    to_pos = chr(nx + ord('A')) + str(3-ny)
                    new_moves = moves + [f"{color},{from_pos},{to_pos}"]
                    
                    queue.append((new_board, new_moves, 'B' if color == 'w' else 'w'))

    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")