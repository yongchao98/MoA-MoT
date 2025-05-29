from collections import deque
import copy

def get_knight_moves(pos):
    moves = [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:  # 4x3 board
            valid.append((new_x, new_y))
    return valid

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at black's initial positions and vice versa
    white_positions = {(0,2), (1,0)}  # Target positions for white
    black_positions = {(0,0), (3,1)}  # Target positions for black
    
    current_white = set()
    current_black = set()
    for i in range(4):
        for j in range(3):
            if board[j][i] == 'w':
                current_white.add((i,j))
            elif board[j][i] == 'B':
                current_black.add((i,j))
    
    return current_white == white_positions and current_black == black_positions

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['B', '.', ' ', '.'],
        [' ', ' ', '.', 'w'],
        ['w', 'B', '.', ' ']
    ]
    
    seen = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        board_state = (board_to_string(board), is_white_turn)
        if board_state in seen:
            continue
        seen.add(board_state)
        
        # Find all pieces that can move
        piece = 'w' if is_white_turn else 'B'
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    # Try all possible moves for this piece
                    for new_x, new_y in get_knight_moves((x, y)):
                        if board[new_y][new_x] == '.':
                            # Make move
                            new_board = copy.deepcopy(board)
                            new_board[y][x] = '.'
                            new_board[new_y][new_x] = piece
                            move = f"{piece},{pos_to_notation((x,y))},{pos_to_notation((new_x,new_y))}"
                            queue.append((new_board, moves + [move], not is_white_turn))

    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")