from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 3 and 0 <= y < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Original positions
    orig_white = {(0,2), (2,1)}  # C3, B1 positions
    orig_black = {(0,0), (2,2)}  # B3, C1 positions
    
    white_positions = set()
    black_positions = set()
    
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                white_positions.add((i,j))
            elif board[i][j] == 'B':
                black_positions.add((i,j))
    
    return white_positions == orig_black and black_positions == orig_white

def convert_to_algebraic(i, j):
    return f"{chr(ord('A')+j)}{3-i}"

def solve_knight_swap():
    # Initial board state (3x4)
    initial_board = [
        ['B', '.', 'w', ' '],
        ['.', ' ', '.', ' '],
        [' ', 'w', 'B', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(board):
            return moves
        
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            i, j = piece_pos
            for new_i, new_j in get_knight_moves(piece_pos):
                if board[new_i][new_j] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[i][j] = '.'
                    new_board[new_i][new_j] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        from_sq = convert_to_algebraic(i, j)
                        to_sq = convert_to_algebraic(new_i, new_j)
                        new_moves = moves + [f"{piece},{from_sq},{to_sq}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

result = solve_knight_swap()
if result is None:
    print("No")
else:
    print(result)