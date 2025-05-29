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
    # Check if white knights are at B1,C1 and black knights are at B3,D1
    white_positions = set(get_piece_positions(board, 'w'))
    black_positions = set(get_piece_positions(board, 'B'))
    target_white = {(1,2), (2,2)}  # B1,C1 in array coordinates
    target_black = {(1,0), (3,2)}  # B3,D1 in array coordinates
    return white_positions == target_white and black_positions == target_black

def print_board(board):
    for row in board:
        print(row)
    print()

def solve_knight_swap():
    # Initial board (array coordinates)
    initial_board = [
        ['.', 'w', '.', '.'],  # row 3
        ['.', '.', '.', '.'],  # row 2
        ['.', 'B', 'B', 'w']   # row 1
    ]
    
    visited = set()
    queue = deque([(initial_board, [], 'w')])
    
    while queue:
        board, moves, color = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        if is_target_state(board):
            return moves
        
        # Get all pieces of current color
        pieces = get_piece_positions(board, color)
        
        # Try all possible moves for each piece
        for px, py in pieces:
            for nx, ny in get_knight_moves((px, py)):
                if board[ny][nx] == '.':  # If target square is empty
                    # Make the move
                    new_board = copy.deepcopy(board)
                    new_board[py][px] = '.'
                    new_board[ny][nx] = color
                    
                    # Convert to chess notation (A1, B2, etc.)
                    from_pos = chr(px + ord('A')) + str(3-py)
                    to_pos = chr(nx + ord('A')) + str(3-ny)
                    new_moves = moves + [f"{color},{from_pos},{to_pos}"]
                    
                    # Add new state to queue
                    queue.append((new_board, new_moves, 'B' if color == 'w' else 'w'))
    
    return ["No"]

# Run the solver
solution = solve_knight_swap()
print(solution)