import json
from collections import deque
from copy import deepcopy

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

def move_box(state, box_id, target, visited_states=None):
    if visited_states is None:
        visited_states = set()
    
    queue = deque([(state, [])])
    visited = {json.dumps(state, sort_keys=True)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state[box_id] == target:
            return path + [current_state]
        
        occupied = set(current_state.values())
        current_pos = current_state[box_id]
        
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                next_state = deepcopy(current_state)
                next_state[box_id] = next_pos
                state_str = json.dumps(next_state, sort_keys=True)
                
                if state_str not in visited and state_str not in visited_states:
                    visited.add(state_str)
                    queue.append((next_state, path + [current_state]))
    
    return None

def solve():
    current_state = deepcopy(initial_state)
    solution = [current_state]
    visited_states = {json.dumps(current_state, sort_keys=True)}
    
    # First move box4 to temporary position to clear path
    path = move_box(current_state, "box4", "C2,4", visited_states)
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
        for state in path:
            visited_states.add(json.dumps(state, sort_keys=True))
    
    # Move box3 to its goal
    path = move_box(current_state, "box3", "C1,2", visited_states)
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
        for state in path:
            visited_states.add(json.dumps(state, sort_keys=True))
    
    # Move box2 to its goal
    path = move_box(current_state, "box2", "C2,2", visited_states)
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
        for state in path:
            visited_states.add(json.dumps(state, sort_keys=True))
    
    # Move box4 to its final goal
    path = move_box(current_state, "box4", "C2,5", visited_states)
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
        for state in path:
            visited_states.add(json.dumps(state, sort_keys=True))
    
    # Finally move box1 to its goal
    path = move_box(current_state, "box1", "C2,8", visited_states)
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
        for state in path:
            visited_states.add(json.dumps(state, sort_keys=True))
    
    # Verify solution
    if current_state == goal_state:
        # Remove consecutive duplicates
        simplified = [solution[0]]
        for state in solution[1:]:
            if state != simplified[-1]:
                simplified.append(state)
        return simplified
    
    return None

# Solve and print result
solution = solve()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")