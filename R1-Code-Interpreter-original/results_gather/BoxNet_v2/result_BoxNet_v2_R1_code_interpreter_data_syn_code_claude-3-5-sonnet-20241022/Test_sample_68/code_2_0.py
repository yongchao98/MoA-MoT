import heapq
from copy import deepcopy
import json

# Problem setup
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

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in current_state.values():
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

def solve():
    def get_state_key(state):
        return tuple(sorted((k, v) for k, v in state.items()))
    
    start = initial_state
    queue = [(0, start)]
    visited = {get_state_key(start): 0}
    parent = {get_state_key(start): None}
    
    while queue:
        cost, current = queue.pop(0)
        
        if current == goal_state:
            # Reconstruct path
            path = []
            current_key = get_state_key(current)
            while current_key is not None:
                state_dict = dict(current_key)
                path.append(state_dict)
                current_key = parent[current_key]
            return path[::-1]
        
        for next_state in get_next_states(current):
            next_key = get_state_key(next_state)
            new_cost = cost + 1
            
            if next_key not in visited or new_cost < visited[next_key]:
                visited[next_key] = new_cost
                parent[next_key] = get_state_key(current)
                queue.append((new_cost, next_state))
        
        queue.sort(key=lambda x: x[0])  # Sort by cost
    
    return None

# Solve and print result
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")