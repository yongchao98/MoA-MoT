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

def get_next_states(state):
    states = []
    for box in state:
        if state[box] == goal_state[box]:
            continue
        curr_pos = state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                states.append(new_state)
    return states

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def solve():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
        
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

# Get solution and print in compact format
solution = solve()
if solution:
    result = []
    for state in solution:
        result.append({
            "box1": state["box1"],
            "box2": state["box2"],
            "box3": state["box3"],
            "box4": state["box4"]
        })
    print(json.dumps(result))
else:
    print("No solution found")