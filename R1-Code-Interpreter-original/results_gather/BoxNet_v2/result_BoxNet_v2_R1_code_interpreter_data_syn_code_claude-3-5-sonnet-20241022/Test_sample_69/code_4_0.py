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

def move_box(state, box, new_pos):
    new_state = dict(state)
    new_state[box] = new_pos
    return new_state

def is_cell_free(pos, state, exclude_box=None):
    return not any(box != exclude_box and state[box] == pos for box in state)

def solve():
    states = [dict(initial_state)]
    current = dict(initial_state)
    
    # Predefined sequence of moves
    moves = [
        # Move box1 to goal (C2,5)
        ("box1", ["C4,5", "C3,5", "C2,5"]),
        # Move box2 to goal (C3,4)
        ("box2", ["C2,3", "C3,3", "C3,4"]),
        # Move box4 to goal (C3,1)
        ("box4", ["C3,5", "C3,4", "C3,3", "C3,2", "C3,1"]),
        # Move box5 to goal (C1,3)
        ("box5", ["C3,4", "C2,4", "C1,4", "C1,3"]),
        # Move box3 to goal (C1,5)
        ("box3", ["C2,3", "C1,3", "C1,4", "C1,5"])
    ]
    
    for box, path in moves:
        if current[box] != goal_state[box]:
            for pos in path:
                if is_cell_free(pos, current, box):
                    current = move_box(current, box, pos)
                    states.append(dict(current))
                    if current[box] == goal_state[box]:
                        break
    
    return states

# Generate and print solution
solution = solve()
print(json.dumps(solution))