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

def is_cell_occupied(state, cell, exclude_box=None):
    for box, pos in state.items():
        if box != exclude_box and pos == cell:
            return True
    return False

def get_path(start, end, current_state, moving_box):
    if start == end:
        return []
    
    visited = {start}
    queue = [(start, [])]
    
    while queue:
        current, path = queue.pop(0)
        for next_cell in adjacency[current]:
            if next_cell not in visited and not is_cell_occupied(current_state, next_cell, moving_box):
                if next_cell == end:
                    return path + [next_cell]
                visited.add(next_cell)
                queue.append((next_cell, path + [next_cell]))
    return None

def solve():
    solution = [dict(initial_state)]
    current_state = dict(initial_state)
    
    # Process boxes in order
    for box in ["box1", "box2", "box3", "box4", "box5"]:
        if current_state[box] == goal_state[box]:
            continue
            
        path = get_path(current_state[box], goal_state[box], current_state, box)
        if path:
            for next_pos in path:
                current_state = dict(current_state)
                current_state[box] = next_pos
                solution.append(dict(current_state))
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))