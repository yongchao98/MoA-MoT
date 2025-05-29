import json
from collections import deque
import copy

initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

def state_to_string(state):
    return str(sorted(state.items()))

def is_valid_move(current_state, box, new_pos):
    # Check if new position is occupied
    for other_box, pos in current_state.items():
        if other_box != box and pos == new_pos:
            return False
    return True

def get_manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_next_states(current_state):
    next_states = []
    # Prioritize boxes that are furthest from their goals
    boxes = list(current_state.keys())
    boxes.sort(key=lambda b: -get_manhattan_distance(current_state[b], goal_state[b]))
    
    for box in boxes:
        if current_state[box] == goal_state[box]:
            continue
        
        current_pos = current_state[box]
        # Sort adjacent positions by distance to goal
        adj_positions = sorted(
            adjacency[current_pos],
            key=lambda pos: get_manhattan_distance(pos, goal_state[box])
        )
        
        for next_pos in adj_positions:
            if is_valid_move(current_state, box, next_pos):
                new_state = copy.deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
        
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

solution = solve()
print(json.dumps(solution))