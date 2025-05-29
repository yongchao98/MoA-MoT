from collections import deque
import copy

def is_valid_knight_move(start, end):
    dx = abs(end[0] - start[0])
    dy = abs(end[1] - start[1])
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_knight_moves(pos):
    moves = [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:
            valid.append((new_x, new_y))
    return valid

def notation_to_pos(notation):
    return (ord(notation[0]) - ord('A'), int(notation[1]) - 1)

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def solve_knight_swap():
    # Initial board state
    board = [
        ['B', '.', ' ', '.'],  # row 3
        [' ', ' ', '.', 'w'],  # row 2
        ['w', 'B', '.', ' ']   # row 1
    ]
    
    initial_white = {(0,2), (3,1)}  # A1, D2
    initial_black = {(0,0), (1,2)}  # A3, B1
    target_white = {(0,0), (1,2)}   # A3, B1
    target_black = {(0,2), (3,1)}   # A1, D2
    
    seen = set()
    queue = deque([(board, [], True, initial_white, initial_black)])
    
    while queue:
        current_board, moves, is_white_turn, white_pos, black_pos = queue.popleft()
        
        if white_pos == target_white and black_pos == target_black:
            return moves
        
        state = (tuple(map(tuple, current_board)), is_white_turn)
        if state in seen:
            continue
        seen.add(state)
        
        current_pieces = white_pos if is_white_turn else black_pos
        piece = 'w' if is_white_turn else 'B'
        
        for start_pos in current_pieces:
            for end_pos in get_knight_moves(start_pos):
                if current_board[end_pos[1]][end_pos[0]] == '.':
                    new_board = copy.deepcopy(current_board)
                    new_board[start_pos[1]][start_pos[0]] = '.'
                    new_board[end_pos[1]][end_pos[0]] = piece
                    
                    new_white = {end_pos if pos == start_pos else pos for pos in white_pos} if is_white_turn else white_pos
                    new_black = black_pos if is_white_turn else {end_pos if pos == start_pos else pos for pos in black_pos}
                    
                    move = f"{piece},{pos_to_notation(start_pos)},{pos_to_notation(end_pos)}"
                    queue.append((new_board, moves + [move], not is_white_turn, new_white, new_black))

    return None

solution = solve_knight_swap()
print(f"<<<{solution}>>>")