from collections import deque
import copy

def get_knight_moves(pos):
    moves = [
        (2, 1), (2, -1), (-2, 1), (-2, -1),
        (1, 2), (1, -2), (-1, 2), (-1, -2)
    ]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:
            valid.append((new_x, new_y))
    return valid

def state_to_string(white_pos, black_pos, white_turn):
    return f"{sorted(white_pos)},{sorted(black_pos)},{white_turn}"

def is_target_reached(white_pos, black_pos):
    target_black = {(0, 1), (1, 1)}
    target_white = {(2, 0), (3, 0)}
    return set(white_pos) == target_white and set(black_pos) == target_black

def find_solution():
    # Initial positions
    white_pos = [(0, 1), (1, 1)]  # A2, B2
    black_pos = [(2, 0), (3, 0)]  # C1, D1
    
    queue = deque([(white_pos, black_pos, True, [])])  # positions, white_turn, moves
    visited = {state_to_string(white_pos, black_pos, True)}
    
    while queue:
        current_white, current_black, white_turn, moves = queue.popleft()
        
        if is_target_reached(current_white, current_black):
            return moves
        
        current_pieces = current_white if white_turn else current_black
        other_pieces = current_black if white_turn else current_white
        
        for i, piece in enumerate(current_pieces):
            for new_pos in get_knight_moves(piece):
                # Check if move is valid (square is empty)
                if new_pos not in current_pieces and new_pos not in other_pieces:
                    # Create new state
                    new_current = list(current_pieces)
                    new_current[i] = new_pos
                    
                    # Create move notation
                    from_sq = f"{chr(piece[0] + 65)}{piece[1] + 1}"
                    to_sq = f"{chr(new_pos[0] + 65)}{new_pos[1] + 1}"
                    move = f"{'w' if white_turn else 'B'},{from_sq},{to_sq}"
                    
                    # Set up next state
                    next_white = new_current if white_turn else current_white
                    next_black = current_black if white_turn else new_current
                    
                    state = state_to_string(next_white, next_black, not white_turn)
                    if state not in visited:
                        visited.add(state)
                        queue.append((next_white, next_black, not white_turn, moves + [move]))
    
    return None

solution = find_solution()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")