from collections import deque

# Define the initial state
initial_state = {
    'player': (7, 5),
    'boxes': {(3, 4), (3, 5), (5, 4), (6, 2), (6, 4)},
    'goals': {(1, 3), (3, 2), (4, 2), (5, 2), (7, 2)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 1), (1, 2), (1, 5), (1, 6),
              (2, 0), (2, 1), (2, 6),
              (3, 0), (3, 1), (3, 6),
              (4, 0), (4, 1), (4, 6),
              (5, 0), (5, 1), (5, 6),
              (6, 0), (6, 6),
              (7, 0), (7, 6),
              (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6)},
    'path': ''
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a move is valid
def is_valid_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    # Check if new player position is a wall
    if new_player_pos in state['walls']:
        return False
    
    # Check if new player position is a box
    if new_player_pos in state['boxes']:
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        # Check if new box position is a wall or another box
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    
    return True

# Function to apply a move
def apply_move(state, move_key, move):
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
        'path': state['path'] + move_key
    }

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# BFS to find the solution
def solve_sokoban(initial_state):
    queue = deque([initial_state])
    visited = set()
    
    while queue:
        current_state = queue.popleft()
        
        # Check if solved
        if is_solved(current_state):
            return current_state['path']
        
        # Mark current state as visited
        state_key = (current_state['player'], frozenset(current_state['boxes']))
        if state_key in visited:
            continue
        visited.add(state_key)
        
        # Try all possible moves
        for move_key, move in moves.items():
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move_key, move)
                queue.append(new_state)
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)