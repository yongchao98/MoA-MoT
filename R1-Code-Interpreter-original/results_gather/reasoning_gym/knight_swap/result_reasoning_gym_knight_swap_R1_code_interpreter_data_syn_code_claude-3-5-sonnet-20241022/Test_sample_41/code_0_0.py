from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(pos[1] + 1)

def state_to_string(state):
    return ','.join([f"{p}" for p in state['white'] + state['black']])

def is_goal_state(state):
    target_white = {(1,2), (0,0)}  # B3, A1
    target_black = {(3,1), (2,0)}  # D2, C1
    return set(state['white']) == target_white and set(state['black']) == target_black

def solve_knight_swap():
    # Initial state
    initial_state = {
        'white': [(3,1), (2,0)],  # D2, C1
        'black': [(1,2), (0,0)],  # B3, A1
        'empty': [(0,2), (1,1), (1,0), (2,1), (3,2), (0,1), (2,2), (3,0)],
        'turn': 'w'
    }
    
    queue = deque([(initial_state, [])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state):
            return moves
        
        pieces = current_state['white'] if current_state['turn'] == 'w' else current_state['black']
        
        for i, piece in enumerate(pieces):
            for new_pos in get_knight_moves(piece):
                if new_pos in current_state['empty']:
                    new_state = copy.deepcopy(current_state)
                    
                    # Update piece position
                    if current_state['turn'] == 'w':
                        new_state['white'][i] = new_pos
                    else:
                        new_state['black'][i] = new_pos
                    
                    # Update empty squares
                    new_state['empty'].remove(new_pos)
                    new_state['empty'].append(piece)
                    
                    # Switch turn
                    new_state['turn'] = 'b' if current_state['turn'] == 'w' else 'w'
                    
                    state_str = state_to_string(new_state)
                    if state_str not in visited:
                        visited.add(state_str)
                        new_moves = moves + [f"{current_state['turn']},{pos_to_notation(piece)},{pos_to_notation(new_pos)}"]
                        queue.append((new_state, new_moves))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")