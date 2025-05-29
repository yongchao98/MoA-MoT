from collections import deque

def is_goal_state(state):
    # Check if all boxes are on goals
    for box in state['boxes']:
        if box not in state['goals']:
            return False
    return True

def get_possible_moves(state):
    # Generate possible moves from the current state
    moves = []
    player = state['player']
    boxes = set(state['boxes'])
    walls = set(state['walls'])
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    for move, (dr, dc) in directions.items():
        new_player = (player[0] + dr, player[1] + dc)
        
        if new_player in walls:
            continue
        
        if new_player in boxes:
            new_box = (new_player[0] + dr, new_player[1] + dc)
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
        
        for move, new_state in get_possible_moves(current_state):
            state_tuple = (new_state['player'], frozenset(new_state['boxes']))
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((new_state, path + move))
    
    return "No solution"

initial_state = {
    'player': (4, 2),
    'boxes': {(1, 3), (3, 2), (4, 3)},
    'goals': {(1, 1), (1, 2), (2, 2)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
              (1, 0), (1, 5), (2, 0), (2, 5), (3, 0), (3, 5),
              (4, 0), (4, 5), (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5)}
}

solution = solve_sokoban(initial_state)
print(solution)