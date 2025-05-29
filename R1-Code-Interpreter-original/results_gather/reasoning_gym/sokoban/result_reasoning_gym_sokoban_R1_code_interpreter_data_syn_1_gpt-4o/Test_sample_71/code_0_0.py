from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (5, 4),
    'boxes': {(3, 3), (5, 3), (5, 5), (6, 5)},
    'goals': {(4, 5), (5, 4), (6, 4), (6, 6)},
    'walls': {(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7),
              (2, 1), (2, 7), (3, 1), (3, 7), (4, 1), (4, 7),
              (5, 1), (5, 7), (6, 1), (6, 7), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6), (7, 7)},
    'empty': {(2, 2), (2, 3), (2, 4), (2, 5), (2, 6),
              (3, 2), (3, 4), (3, 5), (3, 6),
              (4, 2), (4, 3), (4, 4), (4, 6),
              (5, 2), (5, 6),
              (6, 2), (6, 3), (6, 6)}
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
    
    if new_player_pos in state['walls']:
        return False
    
    if new_player_pos in state['boxes']:
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    
    return True

# Function to apply a move
def apply_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    new_boxes = set(state['boxes'])
    if new_player_pos in state['boxes']:
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        new_boxes.remove(new_player_pos)
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls'],
        'empty': state['empty']
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