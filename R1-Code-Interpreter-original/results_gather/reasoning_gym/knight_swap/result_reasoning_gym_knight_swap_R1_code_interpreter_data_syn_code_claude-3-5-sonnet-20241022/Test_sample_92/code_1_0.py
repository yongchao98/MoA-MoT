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
    return [(r, c) for r, c in moves if 0 <= r <= 2 and 0 <= c <= 3]

def board_to_string(board):
    return tuple(tuple(row) for row in board)

def is_target_state(board):
    # White knights should be at original Black positions (B3, D2)
    # Black knights should be at original White positions (B2, C3)
    target_positions = {
        'w': {(0, 1), (1, 3)},  # B3, D2
        'B': {(1, 1), (0, 2)}   # B2, C3
    }
    
    current_positions = {'w': set(), 'B': set()}
    for i in range(3):
        for j in range(4):
            if board[i][j] in ['w', 'B']:
                current_positions[board[i][j]].add((i, j))
    
    return (current_positions['w'] == target_positions['w'] and 
            current_positions['B'] == target_positions['B'])

def get_position_notation(row, col):
    return f"{chr(65 + col)}{3-row}"

def find_solution():
    initial_board = [
        ['.', 'B', 'w', '.'],
        ['.', 'w', '.', 'B'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        current_piece = 'w' if is_white_turn else 'B'
        
        # Find all current player's knights
        knights = []
        for r in range(3):
            for c in range(4):
                if current_board[r][c] == current_piece:
                    knights.append((r, c))
        
        # Try each knight's possible moves
        for kr, kc in knights:
            for new_r, new_c in get_knight_moves((kr, kc)):
                if current_board[new_r][new_c] == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[new_r][new_c] = current_piece
                    new_board[kr][kc] = '.'
                    
                    board_state = board_to_string(new_board)
                    if board_state not in seen:
                        seen.add(board_state)
                        move = f"{current_piece},{get_position_notation(kr, kc)},{get_position_notation(new_r, new_c)}"
                        new_moves = moves + [move]
                        
                        if is_target_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return "No"

result = find_solution()
print(result)