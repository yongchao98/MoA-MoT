from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (3, 7),
    'boxes': {(2, 2), (4, 2), (4, 5), (4, 6), (4, 7)},
    'goals': {(1, 1), (3, 3), (3, 6), (4, 9)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
              (1, 0), (1, 4), (1, 9),
              (2, 0), (2, 8), (2, 9),
              (3, 0), (3, 8), (3, 9),
              (4, 0), (4, 9),
              (5, 0), (5, 9),
              (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 8), (6, 9)},
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
    
    # Check if the new player position is a wall
    if new_player_pos in state['walls']:
        return False
    
    # Check if the new player position is a box
    if new_player_pos in state['boxes']:
        # Calculate the new box position
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        # Check if the new box position is a wall or another box
        if new_box_pos in state['walls'] or new_box_pos in state['boxes']:
            return False
    
    return True

# Function to apply a move
def apply_move(state, move):
    player_pos = state['player']
    new_player_pos = (player_pos[0] + move[0], player_pos[1] + move[1])
    
    new_boxes = set(state['boxes'])
    
    # Check if the new player position is a box
    if new_player_pos in state['boxes']:
        # Calculate the new box position
        new_box_pos = (new_player_pos[0] + move[0], new_player_pos[1] + move[1])
        # Move the box
        new_boxes.remove(new_player_pos)
        new_boxes.add(new_box_pos)
    
    return {
        'player': new_player_pos,
        'boxes': new_boxes,
        'goals': state['goals'],
        'walls': state['walls']
    }

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Breadth-First Search to find the solution
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