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

def get_coordinates(cell):
    row, col = map(int, cell.replace('C', '').split(','))
    return (row, col)

def get_next_move(current_pos, target_pos, occupied_cells):
    current_row, current_col = get_coordinates(current_pos)
    target_row, target_col = get_coordinates(target_pos)
    
    best_move = None
    min_distance = float('inf')
    
    for next_pos in adjacency[current_pos]:
        if next_pos in occupied_cells:
            continue
        next_row, next_col = get_coordinates(next_pos)
        distance = abs(next_row - target_row) + abs(next_col - target_col)
        if distance < min_distance:
            min_distance = distance
            best_move = next_pos
            
    return best_move

def solve():
    solution = [initial_state.copy()]
    current_state = initial_state.copy()
    
    # Process boxes in a specific order
    boxes_order = ["box5", "box2", "box3", "box4", "box1"]
    
    for box in boxes_order:
        while current_state[box] != goal_state[box]:
            occupied_cells = set(current_state.values())
            next_pos = get_next_move(current_state[box], goal_state[box], occupied_cells)
            
            if next_pos:
                new_state = current_state.copy()
                new_state[box] = next_pos
                current_state = new_state
                solution.append(new_state.copy())
    
    print(json.dumps(solution))

solve()