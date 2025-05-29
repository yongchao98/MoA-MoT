import json
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

def find_path(start, end, occupied):
    if start == end:
        return []
    
    queue = [(start, [])]
    visited = {start}
    
    while queue:
        current, path = queue.pop(0)
        
        for next_pos in adjacency[current]:
            if next_pos == end:
                return path + [next_pos]
            
            if next_pos not in visited and next_pos not in occupied:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    
    return None

def solve():
    current_state = deepcopy(initial_state)
    solution = [current_state]
    
    # First move box4 to make space
    occupied = set(current_state.values()) - {current_state["box4"]}
    path = find_path(current_state["box4"], "C2,3", occupied)
    if path:
        for pos in path:
            new_state = deepcopy(current_state)
            new_state["box4"] = pos
            current_state = new_state
            solution.append(current_state)
    
    # Move box3 to its goal
    occupied = set(current_state.values()) - {current_state["box3"]}
    path = find_path(current_state["box3"], goal_state["box3"], occupied)
    if path:
        for pos in path:
            new_state = deepcopy(current_state)
            new_state["box3"] = pos
            current_state = new_state
            solution.append(current_state)
    
    # Move box2 to its goal
    occupied = set(current_state.values()) - {current_state["box2"]}
    path = find_path(current_state["box2"], goal_state["box2"], occupied)
    if path:
        for pos in path:
            new_state = deepcopy(current_state)
            new_state["box2"] = pos
            current_state = new_state
            solution.append(current_state)
    
    # Move box4 to its goal
    occupied = set(current_state.values()) - {current_state["box4"]}
    path = find_path(current_state["box4"], goal_state["box4"], occupied)
    if path:
        for pos in path:
            new_state = deepcopy(current_state)
            new_state["box4"] = pos
            current_state = new_state
            solution.append(current_state)
    
    # Finally move box1 to its goal
    occupied = set(current_state.values()) - {current_state["box1"]}
    path = find_path(current_state["box1"], goal_state["box1"], occupied)
    if path:
        for pos in path:
            new_state = deepcopy(current_state)
            new_state["box1"] = pos
            current_state = new_state
            solution.append(current_state)
    
    # Verify final state matches goal state
    if current_state == goal_state:
        return solution
    return None

# Solve and print result
solution = solve()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")