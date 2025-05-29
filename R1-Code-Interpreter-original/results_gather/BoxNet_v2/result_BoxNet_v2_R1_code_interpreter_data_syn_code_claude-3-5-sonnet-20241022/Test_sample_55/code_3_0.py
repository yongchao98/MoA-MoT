import json
import heapq
from copy import deepcopy

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

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def get_next_states(current_state):
    next_states = []
    occupied = set(current_state.values())
    
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
            
        curr_pos = current_state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied or next_pos == current_state[box]:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve_step_by_step():
    def heuristic(state):
        return sum(manhattan_distance(state[box], goal_state[box]) 
                  for box in state if state[box] != goal_state[box])
    
    visited = set()
    start = (heuristic(initial_state), 0, initial_state, [initial_state])
    queue = [start]
    heapq.heapify(queue)
    
    while queue:
        _, steps, current, path = heapq.heappop(queue)
        state_str = state_to_string(current)
        
        if current == goal_state:
            return path
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state in get_next_states(current):
            if state_to_string(next_state) not in visited:
                h = heuristic(next_state)
                g = steps + 1
                f = g + h
                heapq.heappush(queue, (f, g, next_state, path + [next_state]))
                
        # Limit search depth to prevent infinite loops
        if steps > 50:
            continue
    
    return None

def find_solution():
    path = solve_step_by_step()
    if path:
        # Verify the solution
        for i in range(len(path) - 1):
            current = path[i]
            next_state = path[i + 1]
            # Count differences
            diff_count = sum(1 for box in current if current[box] != next_state[box])
            if diff_count > 1:
                return None
        return path
    return None

# Find and print solution
solution = find_solution()
if solution:
    # Format solution for output
    formatted_solution = []
    for state in solution:
        formatted_solution.append(dict(sorted(state.items())))
    print(json.dumps(formatted_solution))
else:
    print("No solution found")