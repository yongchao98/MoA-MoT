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

def move_box(current_pos, target_pos, current_state):
    # Simple path finding between two positions
    path = []
    current = current_pos
    while current != target_pos:
        if int(current.split(',')[1]) < int(target_pos.split(',')[1]):
            # Move right
            next_pos = f"{current.split(',')[0]},{int(current.split(',')[1]) + 1}"
        elif int(current.split(',')[1]) > int(target_pos.split(',')[1]):
            # Move left
            next_pos = f"{current.split(',')[0]},{int(current.split(',')[1]) - 1}"
        elif current.split(',')[0] != target_pos.split(',')[0]:
            # Move between rows
            next_pos = f"{target_pos.split(',')[0]},{current.split(',')[1]}"
        
        if next_pos in adjacency[current]:
            current = next_pos
            path.append(current)
    
    return path

def solve():
    solution = [dict(initial_state)]
    current_state = dict(initial_state)
    
    # Move boxes one by one to their goals
    move_sequence = ["box1", "box2", "box3", "box4"]
    
    for box in move_sequence:
        if current_state[box] != goal_state[box]:
            path = move_box(current_state[box], goal_state[box], current_state)
            for next_pos in path:
                new_state = dict(current_state)
                new_state[box] = next_pos
                current_state = new_state
                solution.append(new_state)
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))