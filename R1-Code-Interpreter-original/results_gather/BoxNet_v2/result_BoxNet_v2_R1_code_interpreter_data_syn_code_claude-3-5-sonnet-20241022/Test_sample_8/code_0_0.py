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
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def heuristic(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += manhattan_distance(state[box], goal_state[box])
    return total

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    for box in current_state:
        # Skip if box is already at goal
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    start_state = initial_state
    visited = set()
    pq = [(heuristic(start_state), 0, start_state, [start_state])]
    heapq.heapify(pq)
    
    while pq:
        _, cost, current_state, path = heapq.heappop(pq)
        current_string = state_to_string(current_state)
        
        if current_state == goal_state:
            return path
            
        if current_string in visited:
            continue
            
        visited.add(current_string)
        
        for next_state in get_next_states(current_state):
            if state_to_string(next_state) not in visited:
                new_cost = cost + 1
                new_path = path + [next_state]
                priority = new_cost + heuristic(next_state)
                heapq.heappush(pq, (priority, new_cost, next_state, new_path))
    
    return None

# Solve and print result
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")