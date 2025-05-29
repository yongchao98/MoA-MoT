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

def get_state_key(board):
    return tuple(sorted((pos, piece) for pos, piece in board.items()))

def is_goal_state(board, initial_whites, initial_blacks):
    current_whites = {pos for pos, piece in board.items() if piece == 'w'}
    current_blacks = {pos for pos, piece in board.items() if piece == 'B'}
    return current_whites == initial_blacks and current_blacks == initial_whites

def solve_knight_swap():
    # Initial board setup
    board = {
        'A3': 'B', 'B3': 'w', 'C3': '.', 'D3': ' ',
        'A2': 'w', 'B2': '.', 'C2': '.', 'D2': 'B',
        'A1': '.', 'B1': '.', 'C1': '.', 'D1': '.'
    }
    
    initial_whites = {pos for pos, piece in board.items() if piece == 'w'}
    initial_blacks = {pos for pos, piece in board.items() if piece == 'B'}
    
    queue = deque([(board, [], True)])  # board, moves, is_white_turn
    visited = {get_state_key(board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(current_board, initial_whites, initial_blacks):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        positions = [pos for pos, p in current_board.items() if p == piece]
        
        for pos in positions:
            valid_moves = get_knight_moves(pos)
            for new_pos in valid_moves:
                if current_board[new_pos] == '.':
                    new_board = copy.deepcopy(current_board)
                    new_board[new_pos] = piece
                    new_board[pos] = '.'
                    
                    state_key = get_state_key(new_board)
                    if state_key not in visited:
                        visited.add(state_key)
                        new_moves = moves + [f"{piece},{pos},{new_pos}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")