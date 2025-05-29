import json
from collections import deque

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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def get_next_states(state):
    next_states = []
    for box in state:
        if state[box] == goal_state[box]:
            continue
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append((heuristic(new_state), new_state))
    return sorted(next_states)  # Sort by heuristic value

def solve():
    start_state = initial_state
    queue = [(heuristic(start_state), start_state, [start_state])]
    visited = {json.dumps(start_state, sort_keys=True)}
    
    while queue:
        _, current_state, path = queue.pop(0)
        
        if current_state == goal_state:
            return path
        
        for _, next_state in get_next_states(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((heuristic(next_state), next_state, new_path))
        
        queue.sort()  # Sort by heuristic value
    
    return None

# Find and format solution
solution = solve()
if solution:
    # Format each state in the solution
    formatted_solution = []
    for state in solution:
        formatted_state = dict(sorted(state.items()))
        formatted_solution.append(formatted_state)
    
    # Print the complete solution in a single line
    result = json.dumps(formatted_solution)
    print(result)
else:
    print("No solution found")