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
    return [(r,c) for r,c in moves if 0 <= r < 4 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at original black positions and vice versa
    white_positions = set()
    black_positions = set()
    target_white = {(2,1), (1,2)}  # Original black positions
    target_black = {(2,3), (0,2)}  # Original white positions
    
    for i in range(4):
        for j in range(4):
            if board[i][j] == 'w':
                white_positions.add((i,j))
            elif board[i][j] == 'B':
                black_positions.add((i,j))
    
    return white_positions == target_white and black_positions == target_black

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.', 'B', '.', 'w'],
        ['.', '.', 'B', '.'],
        ['.', '.', 'w', '.'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
        
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for i in range(4):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            for next_pos in get_knight_moves(piece_pos):
                if board[next_pos[0]][next_pos[1]] == '.':  # Empty square
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[piece_pos[0]][piece_pos[1]] = '.'
                    new_board[next_pos[0]][next_pos[1]] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Record move in chess notation (converting to A1 format)
                        from_square = f"{chr(piece_pos[1]+65)}{4-piece_pos[0]}"
                        to_square = f"{chr(next_pos[1]+65)}{4-next_pos[0]}"
                        new_moves = moves + [f"{piece},{from_square},{to_square}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
print(solution if solution else "No")