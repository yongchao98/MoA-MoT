from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 4),
    'boxes': {(1, 3), (2, 4), (3, 4)},
    'goals': {(1, 2), (6, 4)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
              (1, 0), (1, 1), (1, 5),
              (2, 0), (2, 1), (2, 2), (2, 5),
              (3, 0), (3, 1), (3, 2), (3, 5),
              (4, 0), (4, 1), (4, 2), (4, 5),
              (5, 0), (5, 1), (5, 2), (5, 5),
              (6, 0), (6, 1), (6, 2), (6, 5),
              (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5)},
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
    player = state['player']
    dx, dy = moves[direction]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in state['walls']:
        return None
    
    new_boxes = set(state['boxes'])
    if new_player in new_boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if new_box in state['walls'] or new_box in new_boxes:
            return None
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    
    return {
        'player': new_player,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls'],
    }

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for direction in moves:
            new_state = move_player(current_state, direction)
            if new_state and (new_state['player'], frozenset(new_state['boxes'])) not in visited:
                visited.add((new_state['player'], frozenset(new_state['boxes'])))
                queue.append((new_state, path + direction))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)