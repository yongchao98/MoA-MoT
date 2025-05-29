from collections import deque
import copy

def get_knight_moves(pos):
    x, y = ord(pos[0]) - ord('A'), int(pos[1]) - 1
    moves = [
        (x-2, y-1), (x-2, y+1), (x-1, y-2), (x-1, y+2),
        (x+1, y-2), (x+1, y+2), (x+2, y-1), (x+2, y+1)
    ]
    valid_moves = []
    for nx, ny in moves:
        if 0 <= nx < 4 and 0 <= ny < 3:
            valid_moves.append(f"{chr(nx + ord('A'))}{ny + 1}")
    return valid_moves

def board_to_string(board):
    return ','.join(f"{pos}:{piece}" for pos, piece in sorted(board.items()) if piece != '.')

def is_target_reached(board):
    return (board.get('A1') == 'w' and board.get('B2') == 'w' and 
            board.get('A2') == 'B' and board.get('D1') == 'B')

def solve_knight_swap():
    # Initial board state
    initial_board = {
        'A1': 'B', 'A2': 'w', 'A3': '.',
        'B1': '.', 'B2': 'B', 'B3': '.',
        'C1': '.', 'C2': '.', 'C3': '.',
        'D1': 'w', 'D2': '.', 'D3': '.'
    }
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = set()
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if (board_state, is_white_turn) in visited:
            continue
        
        visited.add((board_state, is_white_turn))
        
        if is_target_reached(board):
            return moves
        
        # Find all pieces of current color
        current_color = 'w' if is_white_turn else 'B'
        pieces = [pos for pos, piece in board.items() if piece == current_color]
        
        # Try all possible moves for each piece
        for pos in pieces:
            for new_pos in get_knight_moves(pos):
                if board[new_pos] == '.':  # If target square is empty
                    new_board = copy.deepcopy(board)
                    new_board[new_pos] = current_color
                    new_board[pos] = '.'
                    new_moves = moves + [f"{current_color},{pos},{new_pos}"]
                    queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")