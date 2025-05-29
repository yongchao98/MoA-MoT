import json

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

def is_valid_state(state):
    positions = list(state.values())
    return len(positions) == len(set(positions))

def get_path(current_pos, target_pos, occupied_positions):
    if current_pos == target_pos:
        return [current_pos]
        
    visited = {current_pos}
    queue = [(current_pos, [current_pos])]
    
    while queue:
        pos, path = queue.pop(0)
        for next_pos in adjacency[pos]:
            if next_pos not in visited and next_pos not in occupied_positions:
                if next_pos == target_pos:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def solve():
    solution = [dict(initial_state)]
    current_state = dict(initial_state)
    
    # Order boxes by priority (boxes that are blocking others should move first)
    box_order = ["box5", "box2", "box3", "box1", "box4"]
    
    for box in box_order:
        if current_state[box] == goal_state[box]:
            continue
            
        while current_state[box] != goal_state[box]:
            occupied_positions = {pos for b, pos in current_state.items() if b != box}
            path = get_path(current_state[box], goal_state[box], occupied_positions)
            
            if not path or len(path) < 2:
                return None
                
            next_pos = path[1]
            new_state = dict(current_state)
            new_state[box] = next_pos
            
            if is_valid_state(new_state):
                current_state = new_state
                solution.append(dict(current_state))
            else:
                return None
    
    return solution

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")