from collections import deque

# Define the initial state
initial_state = {
    'player': (5, 3),
    'boxes': {(2, 4), (3, 4), (4, 5)},
    'goals': {(1, 5), (4, 3), (5, 5)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 6), (2, 0), (2, 6), (3, 0), (3, 6),
              (4, 0), (4, 6), (5, 0), (5, 6), (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6)},
    'width': 7,
    'height': 7
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

def is_valid_position(pos, state):
    return (0 <= pos[0] < state['height'] and
            0 <= pos[1] < state['width'] and
            pos not in state['walls'])

def move_player(state, direction):
    player = state['player']
    dx, dy = moves[direction]
    new_player = (player[0] + dx, player[1] + dy)
    
    if not is_valid_position(new_player, state):
        return None
    
    new_boxes = set(state['boxes'])
    if new_player in state['boxes']:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if not is_valid_position(new_box, state) or new_box in state['boxes']:
            return None
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    
    return {
        'player': new_player,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls'],
        'width': state['width'],
        'height': state['height']
    }

def bfs(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for direction in moves:
            new_state = move_player(current_state, direction)
            if new_state is None:
                continue
            
            state_id = (new_state['player'], frozenset(new_state['boxes']))
            if state_id not in visited:
                visited.add(state_id)
                queue.append((new_state, path + direction))
    
    return None

solution = bfs(initial_state)
print(solution)