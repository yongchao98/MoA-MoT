import json
from collections import deque

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"], "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def is_valid_move(current_state, next_state):
    differences = sum(1 for k in current_state if current_state[k] != next_state[k])
    if differences != 1:
        return False
    
    # Check for box collision
    positions = list(next_state.values())
    if len(positions) != len(set(positions)):
        return False
    
    # Check if move is to adjacent cell
    for box in current_state:
        if current_state[box] != next_state[box]:
            return next_state[box] in adjacency[current_state[box]]
    return True

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        for next_pos in adjacency[current_state[box]]:
            next_state = current_state.copy()
            next_state[box] = next_pos
            if is_valid_move(current_state, next_state):
                next_states.append(next_state)
    return next_states

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def solve():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    max_depth = 20  # Limit search depth
    
    while queue and len(queue) < 1000:  # Limit queue size
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
            
        if len(path) > max_depth:
            continue
            
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")