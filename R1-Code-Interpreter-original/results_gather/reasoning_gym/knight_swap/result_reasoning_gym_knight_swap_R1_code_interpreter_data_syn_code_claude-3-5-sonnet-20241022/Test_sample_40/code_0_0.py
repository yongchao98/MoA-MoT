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
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at original black positions and vice versa
    white_positions = set()
    black_positions = set()
    target_white = {(2,2), (1,0)}  # C3, A2
    target_black = {(2,1), (0,2)}  # B3, C1
    
    for i in range(3):
        for j in range(3):
            if board[i][j] == 'w':
                white_positions.add((i,j))
            elif board[i][j] == 'B':
                black_positions.add((i,j))
    
    return white_positions == target_white and black_positions == target_black

def get_position_name(pos):
    row, col = pos
    return f"{chr(65+col)}{3-row}"

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.', 'w', 'B'],
        ['B', '.', '.'],
        ['.', '.', 'w']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        for i in range(3):
            for j in range(3):
                if board[i][j] == piece:
                    # Try all possible moves for this piece
                    for new_pos in get_knight_moves((i,j)):
                        ni, nj = new_pos
                        if board[ni][nj] == '.':  # Empty square
                            new_board = copy.deepcopy(board)
                            new_board[i][j] = '.'
                            new_board[ni][nj] = piece
                            
                            board_str = board_to_string(new_board)
                            if board_str not in visited:
                                visited.add(board_str)
                                move = f"{piece},{get_position_name((i,j))},{get_position_name((ni,nj))}"
                                queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")