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
    return ''.join(board[pos] if pos in board else '.' 
                  for row in range(3, 0, -1) 
                  for pos in [chr(col + ord('A')) + str(row) for col in range(4)])

def get_next_states(board, is_white_turn):
    piece = 'w' if is_white_turn else 'B'
    current_positions = [pos for pos, p in board.items() if p == piece]
    empty_squares = [pos for row in range(1, 4) 
                    for pos in [chr(col + ord('A')) + str(row) for col in range(4)]
                    if pos not in board]
    
    states = []
    for pos in current_positions:
        for new_pos in get_knight_moves(pos):
            if new_pos in empty_squares:
                new_board = copy.deepcopy(board)
                new_board[new_pos] = new_board.pop(pos)
                states.append((new_board, new_pos, pos))
    return states

def is_target_reached(board):
    white_positions = {pos for pos, p in board.items() if p == 'w'}
    black_initial = {'A3', 'D2'}
    return white_positions == black_initial

def solve_knight_swap():
    initial_board = {'A3': 'B', 'C3': 'w', 'A2': 'w', 'D2': 'B'}
    queue = deque([(initial_board, [], True)])
    visited = {(board_to_string(initial_board), True)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board) and len(moves) > 0:
            return moves
            
        for new_board, new_pos, old_pos in get_next_states(board, is_white_turn):
            state = (board_to_string(new_board), not is_white_turn)
            if state not in visited:
                visited.add(state)
                piece = 'w' if is_white_turn else 'B'
                move = f"{piece},{old_pos},{new_pos}"
                queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")