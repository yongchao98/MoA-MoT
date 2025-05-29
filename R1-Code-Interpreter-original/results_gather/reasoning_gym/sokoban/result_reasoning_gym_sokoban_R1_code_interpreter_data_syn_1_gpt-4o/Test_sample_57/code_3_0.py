from collections import deque

# Define the initial state
initial_state = {
    'player': (1, 4),
    'boxes': {(2, 3), (3, 4), (4, 4), (5, 3)},
    'goals': {(2, 2), (3, 5), (4, 5), (5, 4)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
              (1, 0), (1, 5),
              (2, 0), (2, 5),
              (3, 0), (3, 5),
              (4, 0), (4, 5),
              (5, 0), (5, 5),
              (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5),
              (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5),
              (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5),
              (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Check if the state is a goal state
def is_goal_state(state):
    return state['boxes'] == state['goals']

# Get the new position after a move
def get_new_position(position, move):
    return (position[0] + move[0], position[1] + move[1])

# Check if a move is valid
def is_valid_move(state, move):
    new_player_pos = get_new_position(state['player'], move)
    if new_player_pos in state['walls']:
        return False
    if new_player_pos in state['boxes']:
        new_box_pos = get_new_position(new_player_pos, move)
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    return True

# Apply a move to the state
def apply_move(state, move):
    new_state = {
        'player': get_new_position(state['player'], move),
        'boxes': set(state['boxes']),
        'goals': state['goals'],
        'walls': state['walls'],
    }
    if new_state['player'] in state['boxes']:
        new_box_pos = get_new_position(new_state['player'], move)
        new_state['boxes'].remove(new_state['player'])
        new_state['boxes'].add(new_box_pos)
    return new_state

# Perform BFS to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for move, delta in moves.items():
            if is_valid_move(current_state, delta):
                new_state = apply_move(current_state, delta)
                state_key = (new_state['player'], frozenset(new_state['boxes']))
                if state_key not in visited:
                    visited.add(state_key)
                    queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)