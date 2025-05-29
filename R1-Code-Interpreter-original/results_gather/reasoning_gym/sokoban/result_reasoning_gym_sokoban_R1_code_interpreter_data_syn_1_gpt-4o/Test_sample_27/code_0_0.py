from collections import deque

def is_goal(state):
    # Check if all boxes are on goals
    return all(pos in state['goals'] for pos in state['boxes'])

def get_next_states(state):
    # Generate all possible next states from the current state
    directions = {'L': (0, -1), 'R': (0, 1), 'U': (-1, 0), 'D': (1, 0)}
    next_states = []
    
    for move, (dr, dc) in directions.items():
        new_player = (state['player'][0] + dr, state['player'][1] + dc)
        
        if state['grid'][new_player[0]][new_player[1]] == '+':
            continue  # Wall, can't move
        
        if new_player in state['boxes']:
            new_box = (new_player[0] + dr, new_player[1] + dc)
            if state['grid'][new_box[0]][new_box[1]] in '+@':
                continue  # Wall or another box, can't push
            
            new_boxes = set(state['boxes'])
            new_boxes.remove(new_player)
            new_boxes.add(new_box)
            next_states.append((move, {'player': new_player, 'boxes': new_boxes, 'goals': state['goals'], 'grid': state['grid']}))
        else:
            next_states.append((move, {'player': new_player, 'boxes': state['boxes'], 'goals': state['goals'], 'grid': state['grid']}))
    
    return next_states

def solve_sokoban(grid, player, boxes, goals):
    initial_state = {'player': player, 'boxes': set(boxes), 'goals': set(goals), 'grid': grid}
    queue = deque([(initial_state, "")])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
        
        state_id = (current_state['player'], frozenset(current_state['boxes']))
        if state_id in visited:
            continue
        visited.add(state_id)
        
        for move, next_state in get_next_states(current_state):
            queue.append((next_state, path + move))
    
    return None

# Define the grid and initial positions
grid = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '@', '-', '-', 'X', '+'],
    ['+', '-', '-', 'X', '@', '-', '+'],
    ['+', '-', 'X', '-', '@', '*', '+'],
    ['+', '-', '-', '-', '$', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

player = (3, 5)
boxes = [(1, 2), (2, 4), (3, 4)]
goals = [(1, 5), (2, 3), (3, 2)]

solution = solve_sokoban(grid, player, boxes, goals)
print(solution)