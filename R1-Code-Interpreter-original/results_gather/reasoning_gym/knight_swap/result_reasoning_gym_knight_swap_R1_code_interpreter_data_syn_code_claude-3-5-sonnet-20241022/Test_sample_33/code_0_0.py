from collections import deque
import copy

def get_knight_moves(pos):
    x, y = ord(pos[0]) - ord('A'), int(pos[1]) - 1
    moves = [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:
            valid.append(chr(new_x + ord('A')) + str(new_y + 1))
    return valid

def board_to_string(board):
    return ','.join(sorted([f"{piece}:{pos}" for pos, piece in board.items() if piece in ['B', 'w']]))

def is_target_state(board):
    return (board.get('B2') == 'B' and board.get('B1') == 'B' and 
            board.get('A1') == 'w' and board.get('D3') == 'w')

def solve_knight_swap():
    # Initial board state
    initial_board = {
        'A1': 'B', 'B1': 'w', 'B2': 'w', 'D3': 'B',
        'A3': '.', 'C1': '.', 'C3': '.', 'D2': '.',
        'A2': '.', 'C2': '.'
    }
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # True means Black's turn
    
    while queue:
        current_board, moves, black_turn = queue.popleft()
        board_state = board_to_string(current_board)
        
        if board_state in visited:
            continue
        
        visited.add(board_state)
        
        if is_target_state(current_board):
            return moves
        
        # Get all pieces of current player
        current_pieces = ['B'] if black_turn else ['w']
        piece_positions = [pos for pos, piece in current_board.items() if piece in current_pieces]
        
        # Try all possible moves for current player's pieces
        for pos in piece_positions:
            piece = current_board[pos]
            for new_pos in get_knight_moves(pos):
                if current_board.get(new_pos) == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[pos] = '.'
                    new_board[new_pos] = piece
                    
                    # Add to queue with move notation
                    move = f"{piece},{pos},{new_pos}"
                    queue.append((new_board, moves + [move], not black_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")