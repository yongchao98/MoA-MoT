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

def is_target_state(board):
    # Check if white knights are at original black positions and vice versa
    white_positions = set()
    black_positions = set()
    target_white = {(1,1), (2,0)}  # B2, C1 positions
    target_black = {(2,2), (1,0)}  # C3, B1 positions
    
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                white_positions.add((i,j))
            elif board[i][j] == 'B':
                black_positions.add((i,j))
    
    return white_positions == target_white and black_positions == target_black

def find_solution():
    # Initial board state
    initial_board = [
        ['.', '.', 'w', '.'],
        ['.', 'B', '.', '.'],
        ['.', 'w', 'B', '.']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
            
        visited.add(board_state)
        
        if is_target_state(board):
            return moves
            
        # Find all pieces of current color
        piece = 'w' if white_turn else 'B'
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for pos in pieces:
            for new_pos in get_knight_moves(pos):
                if board[new_pos[0]][new_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[new_pos[0]][new_pos[1]] = piece
                    new_board[pos[0]][pos[1]] = '.'
                    
                    # Record move
                    from_sq = f"{chr(65+pos[1])}{3-pos[0]}"
                    to_sq = f"{chr(65+new_pos[1])}{3-new_pos[0]}"
                    move = f"{piece},{from_sq},{to_sq}"
                    
                    queue.append((new_board, moves + [move], not white_turn))

    return None

solution = find_solution()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")