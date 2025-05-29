import json
from collections import deque

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Simplified adjacency - only horizontal and vertical moves
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"], "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"], "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C1,1", "C2,2", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"], "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C3,1": ["C2,1", "C3,2", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"], "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,1": ["C3,1", "C4,2"], "C4,2": ["C4,1", "C4,3", "C3,2"], "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"], "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C4,7": ["C4,6", "C4,8", "C3,7"], "C4,8": ["C4,7", "C3,8"]
}

def get_next_states(state):
    states = []
    for box in state:
        if state[box] == goal_state[box]:
            continue
        
        current_pos = state[box]
        occupied = set(state.values())
        occupied.remove(current_pos)
        
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                states.append(new_state)
    return states

def solve_puzzle():
    queue = deque([(initial_state, [initial_state])])
    visited = {json.dumps(initial_state, sort_keys=True)}
    max_steps = 12
    
    while queue:
        current_state, path = queue.popleft()
        
        if len(path) > max_steps:
            continue
        
        if current_state == goal_state:
            return path
        
        for next_state in get_next_states(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
                
                # Early termination if path is too long
                if len(queue) > 1000:
                    return None
    
    return None

# Find solution
solution = solve_puzzle()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")