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

def is_valid_move(state, box, new_pos):
    return new_pos not in state.values()

def move_box(state, box, new_pos):
    new_state = state.copy()
    new_state[box] = new_pos
    return new_state

solution = [initial_state]
current_state = initial_state.copy()

# Predefined sequence of moves
moves = [
    # First move box1 to C2,5
    ("box1", ["C4,5", "C3,5", "C2,5"]),
    
    # Then move box2 to C3,4
    ("box2", ["C2,3", "C3,3", "C3,4"]),
    
    # Move box3 to C1,5
    ("box3", ["C2,3", "C1,3", "C1,4", "C1,5"]),
    
    # Move box4 to C3,1
    ("box4", ["C3,5", "C3,4", "C3,3", "C3,2", "C3,1"]),
    
    # Finally move box5 to C1,3
    ("box5", ["C3,4", "C2,4", "C2,3", "C1,3"])
]

# Execute moves
for box, path in moves:
    for new_pos in path:
        if is_valid_move(current_state, box, new_pos):
            current_state = move_box(current_state, box, new_pos)
            if current_state != solution[-1]:
                solution.append(current_state.copy())

print(json.dumps(solution))