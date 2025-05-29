from queue import PriorityQueue
import json
from copy import deepcopy

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

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def get_next_states(current_state):
    next_states = []
    occupied = set(current_state.values())
    
    # Try moving each box that's not at its goal
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
    # Priority queue with (priority, state_string, path)
    start_state = initial_state
    pq = [(heuristic(start_state), state_to_string(start_state), [start_state])]
    visited = {state_to_string(start_state)}
    
    while pq:
        # Get state with lowest priority
        _, _, path = min(pq)
        pq.remove(min(pq))
        
        current_state = path[-1]
        
        # Check if goal reached
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
        
        # Generate and process next states
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                priority = len(new_path) + heuristic(next_state)
                pq.append((priority, state_str, new_path))
        
        # Keep queue size manageable
        if len(pq) > 5000:
            pq.sort()
            pq = pq[:2500]
    
    return None

# Generate and format solution
solution = solve()
if solution:
    # Format the solution as a list of states
    formatted_solution = []
    for state in solution:
        formatted_solution.append(dict(sorted(state.items())))
    print(json.dumps(formatted_solution))
else:
    print("No solution found")