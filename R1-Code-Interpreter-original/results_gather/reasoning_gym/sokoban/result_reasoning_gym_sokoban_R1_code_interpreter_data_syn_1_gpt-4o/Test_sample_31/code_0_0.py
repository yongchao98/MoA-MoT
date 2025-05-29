from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (7, 4),
    'boxes': {(3, 3), (3, 4), (3, 6), (4, 3), (4, 5), (6, 3)},
    'goals': {(1, 4), (1, 7), (2, 3), (4, 6), (6, 2), (6, 4)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
              (1, 0), (1, 8), (2, 0), (2, 8), (3, 0), (3, 8), (4, 0), (4, 8),
              (5, 0), (5, 8), (6, 0), (6, 8), (7, 0), (7, 8), (8, 0), (8, 1),
              (8, 2), (8, 3), (8, 4), (8, 5), (8, 6), (8, 7), (8, 8)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Check if the puzzle is solved
def is_solved(state):
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
    new_player_pos = get_new_position(state['player'], move)
    new_boxes = set(state['boxes'])
    if new_player_pos in state['boxes']:
        new_box_pos = get_new_position(new_player_pos, move)
        new_boxes.remove(new_player_pos)
        new_boxes.add(new_box_pos)
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls'],
    }

# Perform BFS to find the solution
def bfs(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move, delta in moves.items():
            if is_valid_move(current_state, delta):
                new_state = apply_move(current_state, delta)
                new_state_key = (new_state['player'], frozenset(new_state['boxes']))
                if new_state_key not in visited:
                    visited.add(new_state_key)
                    queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = bfs(initial_state)
print(solution)