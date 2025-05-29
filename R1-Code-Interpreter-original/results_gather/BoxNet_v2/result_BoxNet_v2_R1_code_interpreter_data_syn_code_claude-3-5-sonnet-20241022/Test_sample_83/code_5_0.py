import json
from copy import deepcopy
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

def get_next_states(state):
    states = []
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                states.append(new_state)
    
    return states

def solve_for_box(current_state, target_box):
    if current_state[target_box] == goal_state[target_box]:
        return [current_state]
        
    visited = {str(sorted(current_state.items()))}
    queue = deque([(current_state, [current_state])])
    
    while queue:
        state, path = queue.popleft()
        
        for next_state in get_next_states(state):
            if next_state[target_box] == goal_state[target_box]:
                return path + [next_state]
                
            state_str = str(sorted(next_state.items()))
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
                
        if len(queue) > 1000:
            queue = deque(list(queue)[:500])
    
    return None

def solve():
    current_state = initial_state
    full_path = [current_state]
    
    # Try to solve for each box one at a time
    box_order = ["box1", "box2", "box3", "box4", "box5"]
    
    for box in box_order:
        path = solve_for_box(current_state, box)
        if path is None:
            return None
        
        # Add new states to full path (excluding the first state which is already included)
        full_path.extend(path[1:])
        current_state = path[-1]
    
    return full_path

# Run the solver and print the solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")