from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (5, 3),
    'boxes': {(3, 3), (5, 2), (5, 5), (6, 4)},
    'goals': {(4, 5), (5, 4), (6, 3), (6, 6)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
              (1, 0), (1, 6), (2, 0), (2, 6), (3, 0), (3, 6),
              (4, 0), (4, 6), (5, 0), (5, 6), (6, 0), (6, 6),
              (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6)},
}

# Define possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Function to check if a move is valid
def is_valid_move(state, move):
    px, py = state['player']
    dx, dy = moves[move]
    new_player_pos = (px + dx, py + dy)
    
    # Check if the new player position is a wall
    if new_player_pos in state['walls']:
        return False
    
    # Check if the new player position is a box
    if new_player_pos in state['boxes']:
        # Calculate the new box position
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        # Check if the new box position is a wall or another box
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    
    return True

# Function to apply a move
def apply_move(state, move):
    px, py = state['player']
    dx, dy = moves[move]
    new_player_pos = (px + dx, py + dy)
    
    new_boxes = set(state['boxes'])
    if new_player_pos in state['boxes']:
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        new_boxes.remove(new_player_pos)
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls'],
    }

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                new_state_key = (new_state['player'], frozenset(new_state['boxes']))
                
                if new_state_key not in visited:
                    visited.add(new_state_key)
                    queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)