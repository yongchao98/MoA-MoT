from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (4, 4),
    'boxes': {(2, 1), (3, 4), (6, 4), (8, 2)},
    'goals': {(2, 4), (4, 3), (6, 1), (8, 1)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
              (1, 0), (1, 5),
              (2, 0), (2, 5),
              (3, 0), (3, 5),
              (4, 0), (4, 5),
              (5, 0), (5, 5),
              (6, 0), (6, 5),
              (7, 0), (7, 5),
              (8, 0), (8, 5),
              (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
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
        'walls': state['walls']
    }

# Check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Perform BFS to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
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
solution = solve_sokoban(initial_state)
print(solution)