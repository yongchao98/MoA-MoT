from collections import deque
import copy

def get_valid_moves(board, pos):
    moves = [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:  # Check board boundaries
            if board[2-new_y][new_x] == '.':    # Note: board is indexed from bottom
                valid.append((new_x, new_y))
    return valid

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def get_piece_positions(board):
    white_pos = set()
    black_pos = set()
    for y in range(3):
        for x in range(4):
            if board[y][x] == 'w':
                white_pos.add((x, 2-y))  # Convert to coordinate system
            elif board[y][x] == 'B':
                black_pos.add((x, 2-y))
    return white_pos, black_pos

def solve_knight_swap():
    # Initial board (indexed from bottom to top)
    initial_board = [
        ['w', 'B', '.', ' '],  # row 1
        [' ', ' ', '.', 'w'],  # row 2
        ['B', '.', ' ', '.']   # row 3
    ]
    
    white_pos, black_pos = get_piece_positions(initial_board)
    target_white = {(0,2), (1,0)}  # A3, B1
    target_black = {(0,0), (3,1)}  # A1, D2
    
    seen = set()
    queue = deque([(initial_board, [], True, white_pos, black_pos)])
    
    while queue:
        board, moves, is_white_turn, current_white, current_black = queue.popleft()
        
        # Check if we've reached the target state
        if current_white == target_white and current_black == target_black:
            return moves
            
        state = (tuple(map(tuple, board)), is_white_turn)
        if state in seen:
            continue
        seen.add(state)
        
        # Get current positions to move
        current_pieces = current_white if is_white_turn else current_black
        piece = 'w' if is_white_turn else 'B'
        
        # Try all possible moves for current player
        for start_pos in current_pieces:
            for end_pos in get_valid_moves(board, start_pos):
                # Create new board state
                new_board = copy.deepcopy(board)
                new_board[2-start_pos[1]][start_pos[0]] = '.'
                new_board[2-end_pos[1]][end_pos[0]] = piece
                
                # Update positions
                if is_white_turn:
                    new_white = {end_pos if p == start_pos else p for p in current_white}
                    new_black = current_black
                else:
                    new_white = current_white
                    new_black = {end_pos if p == start_pos else p for p in current_black}
                
                move = f"{piece},{pos_to_notation(start_pos)},{pos_to_notation(end_pos)}"
                queue.append((new_board, moves + [move], not is_white_turn, new_white, new_black))
    
    return None

solution = solve_knight_swap()
print(f"<<<{solution}>>>")