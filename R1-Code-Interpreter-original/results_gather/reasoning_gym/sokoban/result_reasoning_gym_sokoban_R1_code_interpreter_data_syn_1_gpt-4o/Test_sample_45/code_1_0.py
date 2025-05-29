from collections import deque

# Define the initial state
initial_state = {
    'player': (6, 6),
    'boxes': {(2, 2), (5, 3), (5, 4), (5, 5), (6, 5)},
    'goals': {(1, 3), (2, 1), (3, 4), (6, 5)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
              (1, 0), (1, 1), (1, 2), (1, 4), (1, 5), (1, 6), (1, 7),
              (2, 0), (2, 6), (2, 7),
              (3, 0), (3, 1), (3, 2), (3, 5), (3, 6), (3, 7),
              (4, 0), (4, 6), (4, 7),
              (5, 0), (5, 6), (5, 7),
              (6, 0), (6, 7),
              (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6), (7, 7)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Function to check if a position is free (not a wall or a box)
def is_free(position, boxes, walls):
    return position not in walls and position not in boxes

# Function to apply a move
def apply_move(state, move):
    player = state['player']
    boxes = state['boxes']
    dx, dy = moves[move]
    new_player = (player[0] + dx, player[1] + dy)
    
    if new_player in boxes:
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if is_free(new_box, boxes, state['walls']):
            new_boxes = boxes - {new_player} | {new_box}
            return {'player': new_player, 'boxes': new_boxes, 'goals': state['goals'], 'walls': state['walls']}
    elif is_free(new_player, boxes, state['walls']):
        return {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}
    
    return None

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        state, path = queue.popleft()
        
        if is_solved(state):
            return path
        
        for move in moves:
            new_state = apply_move(state, move)
            if new_state and (new_state['player'], frozenset(new_state['boxes'])) not in visited:
                visited.add((new_state['player'], frozenset(new_state['boxes'])))
                queue.append((new_state, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)