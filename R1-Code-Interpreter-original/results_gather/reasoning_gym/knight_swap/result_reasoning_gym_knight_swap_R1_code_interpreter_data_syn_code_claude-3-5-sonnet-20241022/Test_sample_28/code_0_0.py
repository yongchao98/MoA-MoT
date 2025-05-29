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
    return [(r, c) for r, c in moves if 0 <= r < 4 and 0 <= c < 4]

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{pos[0]+1}"

def state_to_string(state):
    return ','.join([f"{p}" for p in state])

def solve_knight_swap():
    # Initial state: (white1_pos, white2_pos, black1_pos, black2_pos)
    initial_state = ((1,0), (0,1), (2,2), (0,3))  # A2,B1,C3,D1
    target_state = ((2,2), (0,3), (1,0), (0,1))   # C3,D1,A2,B1
    
    queue = deque([(initial_state, [], True)])  # state, moves, white_turn
    visited = {state_to_string(initial_state)}
    
    while queue:
        state, moves, white_turn = queue.popleft()
        
        if state == target_state:
            return moves
            
        # Get current positions
        w1, w2, b1, b2 = state
        pieces = [w1, w2] if white_turn else [b1, b2]
        
        # Try moving each piece
        for i, piece in enumerate(pieces):
            for new_pos in get_knight_moves(piece):
                # Check if position is empty
                if new_pos not in (w1, w2, b1, b2):
                    # Create new state
                    new_state = list(state)
                    new_state[i if white_turn else i+2] = new_pos
                    new_state = tuple(new_state)
                    
                    state_str = state_to_string(new_state)
                    if state_str not in visited:
                        visited.add(state_str)
                        piece_notation = pos_to_notation(piece)
                        new_pos_notation = pos_to_notation(new_pos)
                        move = f"{'w' if white_turn else 'B'},{piece_notation},{new_pos_notation}"
                        queue.append((new_state, moves + [move], not white_turn))
    
    return None

# Solve and print result
result = solve_knight_swap()
if result:
    print(result)
else:
    print("No")