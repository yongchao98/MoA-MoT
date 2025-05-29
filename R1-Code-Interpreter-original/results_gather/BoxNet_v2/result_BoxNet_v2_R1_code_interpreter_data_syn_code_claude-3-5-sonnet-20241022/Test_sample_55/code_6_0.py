import json
from copy import deepcopy

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

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def is_valid_move(current_state, next_state):
    differences = 0
    for box in current_state:
        if current_state[box] != next_state[box]:
            differences += 1
            if differences > 1:
                return False
            if next_state[box] not in adjacency[current_state[box]]:
                return False
    return differences == 1

def solve_for_box(current_state, box_id, target_pos, visited, path, depth=0):
    if depth > 20:  # Limit depth for single box movement
        return None
    
    if current_state[box_id] == target_pos:
        return path
    
    state_str = str(sorted(current_state.items()))
    if state_str in visited:
        return None
    
    visited.add(state_str)
    
    current_pos = current_state[box_id]
    occupied = set(current_state.values())
    
    # Try all possible moves for the current box
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied or next_pos == current_pos:
            next_state = deepcopy(current_state)
            next_state[box_id] = next_pos
            
            if is_valid_move(current_state, next_state):
                result = solve_for_box(next_state, box_id, target_pos, visited, path + [next_state], depth + 1)
                if result:
                    return result
    
    visited.remove(state_str)
    return None

def solve_complete():
    current_state = deepcopy(initial_state)
    solution = [current_state]
    
    # Order boxes by priority (closest to their goals first)
    boxes = list(current_state.keys())
    boxes.sort(key=lambda box: manhattan_distance(current_state[box], goal_state[box]))
    
    for box in boxes:
        if current_state[box] != goal_state[box]:
            visited = set()
            path = solve_for_box(current_state, box, goal_state[box], visited, [])
            
            if path:
                solution.extend(path[1:])  # Skip the first state as it's already included
                current_state = path[-1]
            else:
                return None
    
    return solution

# Find and verify solution
solution = solve_complete()

if solution:
    # Verify solution
    valid = True
    for i in range(len(solution) - 1):
        if not is_valid_move(solution[i], solution[i+1]):
            valid = False
            break
    
    if valid and solution[-1] == goal_state:
        formatted_solution = []
        for state in solution:
            formatted_solution.append(dict(sorted(state.items())))
        print(json.dumps(formatted_solution))
    else:
        print("No valid solution found")
else:
    print("No solution found")