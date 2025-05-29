from collections import deque

# Define the initial state
initial_state = {
    'player': (3, 4),
    'boxes': {(1, 3), (2, 4), (3, 2), (3, 3)},
    'goals': {(1, 5), (2, 2), (3, 6), (4, 2)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
              (1, 0), (1, 8), (2, 0), (2, 8), (3, 0), (3, 8), (4, 0), (4, 8),
              (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8)},
    'width': 9,
    'height': 6
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

def move_player(state, direction):
    dx, dy = moves[direction]
    px, py = state['player']
    new_player_pos = (px + dx, py + dy)
    
    if new_player_pos in state['walls']:
        return None
    
    new_boxes = set(state['boxes'])
    if new_player_pos in new_boxes:
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        if new_box_pos in state['walls'] or new_box_pos in new_boxes:
            return None
        new_boxes.remove(new_player_pos)
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
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