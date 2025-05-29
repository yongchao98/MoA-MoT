from collections import deque

# Define the initial state
initial_state = {
    'player': (2, 6),
    'boxes': {(3, 4), (4, 1), (4, 6)},
    'goals': {(1, 4), (2, 1), (3, 6)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
              (1, 0), (1, 8),
              (2, 0), (2, 4), (2, 8),
              (3, 0), (3, 8),
              (4, 0), (4, 8),
              (5, 0), (5, 8),
              (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 8)},
}

# Define possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Function to check if a move is valid
def is_valid_move(state, move):
    px, py = state['player']
    dx, dy = moves[move]
    new_pos = (px + dx, py + dy)
    
    if new_pos in state['walls']:
        return False
    
    if new_pos in state['boxes']:
        box_new_pos = (new_pos[0] + dx, new_pos[1] + dy)
        if box_new_pos in state['walls'] or box_new_pos in state['boxes']:
            return False
    
    return True

# Function to apply a move
def apply_move(state, move):
    px, py = state['player']
    dx, dy = moves[move]
    new_pos = (px + dx, py + dy)
    
    new_boxes = set(state['boxes'])
    if new_pos in new_boxes:
        new_boxes.remove(new_pos)
        new_boxes.add((new_pos[0] + dx, new_pos[1] + dy))
    
    return {
        'player': new_pos,
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
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move in moves:
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                state_key = (new_state['player'], frozenset(new_state['boxes']))
                
                if state_key not in visited:
                    visited.add(state_key)
                    queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)