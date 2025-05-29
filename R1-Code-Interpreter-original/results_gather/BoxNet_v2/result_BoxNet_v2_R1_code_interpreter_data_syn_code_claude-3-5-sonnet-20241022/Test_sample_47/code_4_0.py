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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def solve_single_box(state, box_id, target_pos, fixed_positions):
    if state[box_id] == target_pos:
        return [state]
    
    visited = {json.dumps(state, sort_keys=True)}
    queue = [(0, [state])]
    
    while queue:
        _, path = heapq.heappop(queue)
        current = path[-1]
        
        if current[box_id] == target_pos:
            return path
            
        curr_pos = current[box_id]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in fixed_positions:
                new_state = current.copy()
                new_state[box_id] = next_pos
                state_str = json.dumps(new_state, sort_keys=True)
                
                if state_str not in visited:
                    visited.add(state_str)
                    new_path = path + [new_state]
                    priority = len(new_path)
                    heapq.heappush(queue, (priority, new_path))
    return None

def solve():
    current_state = initial_state.copy()
    final_path = [current_state]
    fixed_positions = set()
    
    # Order boxes by priority (minimize conflicts)
    box_order = ["box3", "box4", "box2", "box1"]
    
    for box in box_order:
        target = goal_state[box]
        if current_state[box] != target:
            # Remove fixed positions that will be freed up
            temp_fixed = {pos for b, pos in current_state.items() 
                         if b != box and pos not in fixed_positions}
            
            path = solve_single_box(current_state, box, target, temp_fixed)
            if path is None:
                return None
            
            # Update current state and path
            current_state = path[-1]
            final_path.extend(path[1:])
            fixed_positions.add(target)
    
    return final_path

# Solve and print result
solution = solve()
if solution:
    # Remove duplicate states and ensure proper formatting
    unique_solution = [dict(t) for t in {tuple(sorted(d.items())) for d in solution}]
    unique_solution.sort(key=lambda x: sum(1 for b in x if x[b] != goal_state[b]))
    print("<<<" + json.dumps(unique_solution) + ">>>")
else:
    print("No solution found")