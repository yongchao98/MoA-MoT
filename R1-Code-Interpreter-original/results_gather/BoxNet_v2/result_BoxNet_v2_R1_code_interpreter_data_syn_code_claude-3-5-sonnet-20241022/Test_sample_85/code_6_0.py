import json
from copy import deepcopy
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

def find_path_for_box(start_state, box_id, target_pos, max_steps=20):
    queue = deque([(start_state, [])])
    visited = {json.dumps(start_state, sort_keys=True)}
    
    while queue and len(visited) < 1000:
        current_state, path = queue.popleft()
        
        if current_state[box_id] == target_pos:
            return path + [current_state]
            
        if len(path) >= max_steps:
            continue
            
        occupied = set(current_state.values())
        current_pos = current_state[box_id]
        
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                next_state = deepcopy(current_state)
                next_state[box_id] = next_pos
                state_str = json.dumps(next_state, sort_keys=True)
                
                if state_str not in visited:
                    visited.add(state_str)
                    queue.append((next_state, path + [current_state]))
    
    return None

def solve():
    # Define sequence of moves
    move_sequence = [
        # First clear the path for box3
        ("box4", "C2,3"),  # Move box4 out of the way
        ("box3", "C1,2"),  # Move box3 to its goal
        ("box4", "C2,4"),  # Move box4 towards its goal
        ("box2", "C2,2"),  # Move box2 to its goal
        ("box4", "C2,5"),  # Move box4 to its goal
        ("box1", "C2,8"),  # Move box1 to its goal
    ]
    
    current_state = deepcopy(initial_state)
    full_solution = [current_state]
    
    for box_id, target in move_sequence:
        path = find_path_for_box(current_state, box_id, target)
        if path is None:
            return None
        
        # Add all intermediate states except the first one (to avoid duplication)
        if path:
            full_solution.extend(path[1:])
            current_state = path[-1]
    
    # Verify final state matches goal state
    final_state = full_solution[-1]
    boxes_at_goal = all(
        final_state[box] == goal_state[box] 
        for box in goal_state
    )
    
    if boxes_at_goal:
        # Remove duplicate consecutive states
        simplified = [full_solution[0]]
        for state in full_solution[1:]:
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