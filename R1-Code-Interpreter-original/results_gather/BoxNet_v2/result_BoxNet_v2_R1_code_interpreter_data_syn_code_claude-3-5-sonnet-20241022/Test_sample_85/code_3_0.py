import heapq
from copy import deepcopy
import json

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

def move_single_box(state, box_to_move, target):
    current_path = [state]
    current_state = deepcopy(state)
    visited = {json.dumps(state, sort_keys=True)}
    
    queue = [(0, current_state, current_path)]
    
    while queue:
        _, curr_state, path = queue.pop(0)
        
        if curr_state[box_to_move] == target:
            return path
            
        curr_pos = curr_state[box_to_move]
        occupied = set(curr_state.values())
        
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied:
                next_state = deepcopy(curr_state)
                next_state[box_to_move] = next_pos
                state_str = json.dumps(next_state, sort_keys=True)
                
                if state_str not in visited:
                    visited.add(state_str)
                    new_path = path + [next_state]
                    queue.append((len(new_path), next_state, new_path))
        
        queue.sort(key=lambda x: x[0])
    
    return None

def solve():
    # Define the order of boxes to move
    box_order = [
        ("box4", "C2,1", "C2,2"),  # Move box4 out of the way first
        ("box3", "C2,2", "C1,2"),  # Move box3 to its goal
        ("box2", "C2,8", "C2,2"),  # Move box2 to its goal
        ("box4", "C2,2", "C2,5"),  # Move box4 to its goal
        ("box1", "C1,4", "C2,8")   # Move box1 to its goal
    ]
    
    current_state = deepcopy(initial_state)
    full_path = [current_state]
    
    for box, start, end in box_order:
        if current_state[box] != end:
            path = move_single_box(current_state, box, end)
            if path:
                full_path.extend(path[1:])  # Skip the first state as it's already included
                current_state = path[-1]
            else:
                return None
    
    return full_path

# Solve and print result
solution = solve()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")