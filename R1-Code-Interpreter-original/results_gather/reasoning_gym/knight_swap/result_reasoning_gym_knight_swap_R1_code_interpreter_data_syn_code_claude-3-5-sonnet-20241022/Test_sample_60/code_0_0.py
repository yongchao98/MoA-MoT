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
    # Check if white knights are at original black positions and vice versa
    black_positions = {(2,2), (0,1)}  # Original black positions
    white_positions = {(0,2), (0,3)}  # Original white positions
    
    current_whites = set()
    current_blacks = set()
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                current_whites.add((i,j))
            elif board[i][j] == 'B':
                current_blacks.add((i,j))
    
    return current_whites == black_positions and current_blacks == white_positions

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.', '.', 'B', '.'],
        ['.', '.', '.', '.'],
        ['.', 'B', 'w', 'w']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(current_board):
            return moves
        
        # Find all pieces of current color
        pieces = []
        for i in range(3):
            for j in range(4):
                if current_board[i][j] == ('w' if is_white_turn else 'B'):
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            for new_pos in get_knight_moves(piece_pos):
                if current_board[new_pos[0]][new_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[new_pos[0]][new_pos[1]] = new_board[piece_pos[0]][piece_pos[1]]
                    new_board[piece_pos[0]][piece_pos[1]] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{'w' if is_white_turn else 'B'},{chr(65+piece_pos[1])}{3-piece_pos[0]},{chr(65+new_pos[1])}{3-new_pos[0]}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")