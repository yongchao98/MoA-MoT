import json
from collections import deque

initial_state = {
    "box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"
}
goal_state = {
    "box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C4,2"
}
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"], "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C3,6"]
}

def get_next_states(state):
    states = []
    occupied = set(state.values())
    for box in state:
        if state[box] != goal_state[box]:  # Only move boxes not at goal
            for adj in adjacency[state[box]]:
                if adj not in occupied:
                    new_state = state.copy()
                    new_state[box] = adj
                    states.append(new_state)
    return states

def solve():
    queue = deque([(initial_state, [initial_state])])
    seen = {json.dumps(initial_state, sort_keys=True)}
    
    while queue:
        state, path = queue.popleft()
        if state == goal_state:
            return path
            
        for next_state in get_next_states(state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in seen:
                seen.add(state_str)
                queue.append((next_state, path + [next_state]))
    return None

solution = solve()
if solution:
    print(json.dumps(solution))