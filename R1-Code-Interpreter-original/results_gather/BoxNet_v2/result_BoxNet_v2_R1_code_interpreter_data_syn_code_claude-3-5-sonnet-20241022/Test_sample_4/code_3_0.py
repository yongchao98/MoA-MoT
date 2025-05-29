import copy

initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"], "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"], "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"], "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"], "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C3,6"]
}

def get_next_move(current_pos, target_pos, occupied):
    if current_pos == target_pos:
        return None
    
    visited = {current_pos}
    queue = [(current_pos, [])]
    
    while queue:
        pos, path = queue.pop(0)
        for next_pos in adjacency[pos]:
            if next_pos not in visited and next_pos not in occupied:
                new_path = path + [next_pos]
                if next_pos == target_pos:
                    return new_path[0]
                visited.add(next_pos)
                queue.append((next_pos, new_path))
    return None

def solve():
    solution = [initial_state]
    current = copy.deepcopy(initial_state)
    
    # Order of boxes to move
    box_order = ["box1", "box2", "box3", "box4"]
    completed_boxes = set()
    
    while len(completed_boxes) < 4:
        for box in box_order:
            if box in completed_boxes:
                continue
                
            if current[box] == goal_state[box]:
                completed_boxes.add(box)
                continue
                
            occupied = {pos for b, pos in current.items() if b != box}
            next_pos = get_next_move(current[box], goal_state[box], occupied)
            
            if next_pos:
                new_state = copy.deepcopy(current)
                new_state[box] = next_pos
                current = new_state
                solution.append(new_state)
                break
    
    return solution

# Generate and print solution
result = solve()
print(str(result).replace("'", '"'))