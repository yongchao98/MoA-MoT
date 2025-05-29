from collections import deque

def is_goal_state(state):
    # Check if all boxes are on goals
    for box in state['boxes']:
        if box not in state['goals']:
            return False
    return True

def get_possible_moves(state):
    # Generate possible moves for the player
    moves = []
    player = state['player']
    boxes = set(state['boxes'])
    walls = set(state['walls'])
    
    # Define move directions and their corresponding deltas
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for move, (dx, dy) in directions.items():
        new_player = (player[0] + dx, player[1] + dy)
        
        if new_player in walls:
            continue
        
        if new_player in boxes:
            new_box = (new_player[0] + dx, new_player[1] + dy)
            if new_box in walls or new_box in boxes:
                continue
            new_boxes = set(boxes)
            new_boxes.remove(new_player)
            new_boxes.add(new_box)
            moves.append((move, {'player': new_player, 'boxes': new_boxes, 'goals': state['goals'], 'walls': state['walls']}))
        else:
            moves.append((move, {'player': new_player, 'boxes': boxes, 'goals': state['goals'], 'walls': state['walls']}))
    
    return moves

def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
        
        state_key = (current_state['player'], frozenset(current_state['boxes']))
        if state_key in visited:
            continue
        visited.add(state_key)
        
        for move, new_state in get_possible_moves(current_state):
            queue.append((new_state, path + move))
    
    return "No solution"

# Define the initial state
initial_state = {
    'player': (1, 4),
    'boxes': {(1, 3), (1, 5), (2, 4)},
    'goals': {(3, 3), (5, 3), (3, 5)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
              (1, 0), (1, 1), (1, 2), (1, 6), (1, 7),
              (2, 0), (2, 6), (2, 7),
              (3, 0), (3, 6), (3, 7),
              (4, 0), (4, 1), (4, 6), (4, 7),
              (5, 0), (5, 1), (5, 6), (5, 7),
              (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7)}
}

solution = solve_sokoban(initial_state)
print(solution)