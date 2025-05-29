from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (3, 4),
    'boxes': {(2, 4), (6, 2), (6, 3)},
    'goals': {(2, 2), (4, 3), (4, 4)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 6), (2, 0), (2, 6), (3, 0), (3, 6), (4, 0), (4, 6),
              (5, 0), (5, 6), (6, 0), (6, 6), (7, 0), (7, 6), (8, 0), (8, 1),
              (8, 2), (8, 3), (8, 4), (8, 5), (8, 6)},
    'empty': {(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 1), (2, 3), (2, 5),
              (3, 1), (3, 2), (3, 3), (3, 5), (4, 1), (4, 2), (4, 5), (5, 1),
              (5, 2), (5, 3), (5, 4), (5, 5), (6, 1), (6, 4), (6, 5), (7, 1),
              (7, 2), (7, 3), (7, 4), (7, 5)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Check if the move is valid
def is_valid_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    if new_player_pos in state['walls']:
        return False
    
    if new_player_pos in state['boxes']:
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    
    return True

# Apply the move to the state
def apply_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    new_boxes = set(state['boxes'])
    if new_player_pos in new_boxes:
        new_boxes.remove(new_player_pos)
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls'],
        'empty': state['empty'],
    }

# Check if the state is a goal state
def is_goal_state(state):
    return state['boxes'] == state['goals']

# Perform BFS to find the solution
def bfs(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        for move_name, move in moves.items():
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                state_key = (new_state['player'], frozenset(new_state['boxes']))
                
                if state_key not in visited:
                    visited.add(state_key)
                    queue.append((new_state, path + move_name))
    
    return None

# Find the solution
solution = bfs(initial_state)
print(solution)