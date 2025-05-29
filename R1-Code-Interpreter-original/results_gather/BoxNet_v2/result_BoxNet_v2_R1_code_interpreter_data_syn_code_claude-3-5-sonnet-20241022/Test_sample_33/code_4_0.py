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

def move_box_to_goal(state, box_id, path):
    current = state[box_id]
    target = goal_state[box_id]
    
    if current == target:
        return path
    
    current_state = deepcopy(state)
    while current_state[box_id] != target:
        occupied = set(current_state.values())
        current_pos = current_state[box_id]
        
        # Find best next move towards goal
        best_move = None
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                best_move = next_pos
                if next_pos == target:
                    break
        
        if best_move is None:
            return None
            
        new_state = deepcopy(current_state)
        new_state[box_id] = best_move
        path.append(new_state)
        current_state = new_state
    
    return path

def solve():
    # Move boxes one by one to their goals
    current_state = initial_state
    path = [current_state]
    
    # First move box1 to its goal
    path = move_box_to_goal(current_state, "box1", path)
    if not path:
        return None
    current_state = path[-1]
    
    # Then move box2 to its goal
    path = move_box_to_goal(current_state, "box2", path)
    if not path:
        return None
    current_state = path[-1]
    
    # Then move box3 to its goal
    path = move_box_to_goal(current_state, "box3", path)
    if not path:
        return None
    current_state = path[-1]
    
    # Finally move box4 to its goal
    path = move_box_to_goal(current_state, "box4", path)
    if not path:
        return None
    
    return path

solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")