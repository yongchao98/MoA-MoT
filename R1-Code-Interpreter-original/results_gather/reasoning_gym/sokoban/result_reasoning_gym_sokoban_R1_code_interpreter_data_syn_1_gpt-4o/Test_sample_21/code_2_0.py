from collections import deque

# Define the initial state
initial_state = {
    'player': (1, 1),
    'boxes': {(1, 2), (2, 3), (3, 3), (4, 3)},
    'goals': {(1, 3), (3, 4), (4, 5)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
              (1, 0), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9),
              (2, 0), (2, 4), (2, 5), (2, 8), (2, 9),
              (3, 0), (3, 8), (3, 9),
              (4, 0), (4, 8), (4, 9),
              (5, 0), (5, 8), (5, 9),
              (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 8), (6, 9)},
    'empty': {(1, 5), (1, 6), (2, 2), (2, 6), (3, 1), (3, 2), (3, 5), (3, 6), (3, 7),
              (4, 1), (4, 2), (4, 4), (4, 6), (4, 7), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7)}
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(state):
    return state['boxes'] == state['goals']

def get_new_position(position, move):
    return (position[0] + move[0], position[1] + move[1])

def is_valid_move(state, move):
    new_player_pos = get_new_position(state['player'], move)
    if new_player_pos in state['walls']:
        return False
    if new_player_pos in state['boxes']:
        new_box_pos = get_new_position(new_player_pos, move)
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    return True

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
        'empty': state['empty']
    }

def bfs(initial_state):
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
                state_id = (new_state['player'], frozenset(new_state['boxes']))
                if state_id not in visited:
                    visited.add(state_id)
                    queue.append((new_state, path + move))
    
    return None

solution = bfs(initial_state)
print(solution)