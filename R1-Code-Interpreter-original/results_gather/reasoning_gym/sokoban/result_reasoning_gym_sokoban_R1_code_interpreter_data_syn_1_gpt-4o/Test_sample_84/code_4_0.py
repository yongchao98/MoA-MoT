from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = {
    'player': (7, 5),
    'boxes': {(1, 4), (2, 5), (3, 3), (4, 5), (5, 5), (7, 4)},
    'goals': {(1, 2), (2, 3), (2, 8), (4, 4), (6, 8), (7, 7), (8, 8)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
              (1, 0), (1, 1), (1, 8), (1, 9),
              (2, 0), (2, 1), (2, 9),
              (3, 0), (3, 9),
              (4, 0), (4, 9),
              (5, 0), (5, 9),
              (6, 0), (6, 9),
              (7, 0), (7, 9),
              (8, 0), (8, 9),
              (9, 0), (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 6), (9, 7), (9, 8), (9, 9)},
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

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Improved heuristic function for A* search
def heuristic(state):
    # Sum of Manhattan distances from each box to the nearest goal
    return sum(min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals']) for box in state['boxes'])

# Function to detect deadlocks
def is_deadlock(state):
    # Simple deadlock detection: check if any box is in a corner not on a goal
    for box in state['boxes']:
        if box not in state['goals']:
            if ((box[0] - 1, box[1]) in state['walls'] or (box[0] + 1, box[1]) in state['walls']) and \
               ((box[0], box[1] - 1) in state['walls'] or (box[0], box[1] + 1) in state['walls']):
                return True
    return False

# A* Search to find the solution
def solve_sokoban(initial_state):
    open_set = []
    heappush(open_set, (0, (initial_state['player'], frozenset(initial_state['boxes'])), ""))
    visited = set()
    
    while open_set:
        _, state_key, path = heappop(open_set)
        current_state = {
            'player': state_key[0],
            'boxes': set(state_key[1]),
            'goals': initial_state['goals'],
            'walls': initial_state['walls']
        }
        
        if is_solved(current_state):
            return path
        
        if is_deadlock(current_state):
            continue
        
        for move_name, move in moves.items():
            if is_valid_move(current_state, move):
                new_state = apply_move(current_state, move)
                new_state_key = (new_state['player'], frozenset(new_state['boxes']))
                
                if new_state_key not in visited:
                    visited.add(new_state_key)
                    cost = len(path) + 1 + heuristic(new_state)
                    heappush(open_set, (cost, new_state_key, path + move_name))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)