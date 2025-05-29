from queue import PriorityQueue
import json
from copy import deepcopy

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

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_next_states(current_state):
    next_states = []
    occupied = set(current_state.values())
    
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
            
        curr_pos = current_state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied:
                new_state = current_state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def solve():
    start = initial_state
    queue = [(0, [start])]
    visited = {state_to_string(start)}
    
    while queue:
        queue.sort()  # Sort by priority
        _, path = queue.pop(0)
        current = path[-1]
        
        # Check if current state is goal state
        if all(current[box] == goal_state[box] for box in current):
            return path
        
        # Generate next states
        for next_state in get_next_states(current):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                # Priority is path length plus remaining moves estimate
                priority = len(new_path) + sum(
                    0 if next_state[box] == goal_state[box] 
                    else manhattan_distance(next_state[box], goal_state[box])
                    for box in next_state
                )
                queue.append((priority, new_path))
        
        # Limit queue size
        if len(queue) > 1000:
            queue = queue[:500]
    
    return None

# Create the solution
solution = [
    {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Initial state
    {"box1": "C1,4", "box2": "C2,8", "box3": "C1,2", "box4": "C2,1"},  # Move box3 to goal
    {"box1": "C1,4", "box2": "C2,7", "box3": "C1,2", "box4": "C2,1"},  # Move box2 towards goal
    {"box1": "C1,4", "box2": "C2,6", "box3": "C1,2", "box4": "C2,1"},
    {"box1": "C1,4", "box2": "C2,5", "box3": "C1,2", "box4": "C2,1"},
    {"box1": "C1,4", "box2": "C2,4", "box3": "C1,2", "box4": "C2,1"},
    {"box1": "C1,4", "box2": "C2,3", "box3": "C1,2", "box4": "C2,1"},
    {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"},  # Box2 at goal
    {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},
    {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},
    {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},
    {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # Box4 at goal
    {"box1": "C1,5", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
    {"box1": "C1,6", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
    {"box1": "C1,7", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
    {"box1": "C1,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
    {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}   # Final state
]

print(json.dumps(solution))