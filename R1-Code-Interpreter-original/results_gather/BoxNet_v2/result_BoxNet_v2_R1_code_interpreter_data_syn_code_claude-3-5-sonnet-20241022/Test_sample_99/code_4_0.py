import heapq
from copy import deepcopy
import json

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

def get_next_states(state):
    next_states = []
    occupied = set(state.values())
    
    # First try to move boxes that are not in their goal positions
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def solve():
    visited = set()
    start = (0, 0, [initial_state], initial_state)
    queue = [start]
    heapq.heapify(queue)
    
    while queue:
        _, steps, path, current = heapq.heappop(queue)
        state_str = state_to_string(current)
        
        if current == goal_state:
            return path
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state in get_next_states(current):
            if state_to_string(next_state) not in visited:
                h_score = sum(manhattan_distance(next_state[box], goal_state[box]) 
                            for box in next_state)
                priority = steps + 1 + h_score
                heapq.heappush(queue, (priority, steps + 1, path + [next_state], next_state))
    
    return None

# Solve and print result
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")