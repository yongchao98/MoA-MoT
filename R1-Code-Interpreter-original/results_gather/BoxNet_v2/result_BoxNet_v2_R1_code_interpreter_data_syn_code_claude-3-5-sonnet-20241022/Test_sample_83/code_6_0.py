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

def is_valid_move(state, box, new_pos):
    # Check if the new position is occupied by another box
    return new_pos not in state.values()

def move_box(state, box, new_pos):
    new_state = state.copy()
    new_state[box] = new_pos
    return new_state

def generate_solution():
    solution = [initial_state]
    current_state = initial_state.copy()
    
    # Predefined sequence of moves
    moves = [
        # Move box1 to goal (C2,5)
        ("box1", ["C5,4", "C4,4", "C3,4", "C2,4", "C2,5"]),
        
        # Move box2 to goal (C3,4)
        ("box2", ["C1,4", "C2,4", "C3,4"]),
        
        # Move box3 to goal (C1,5)
        ("box3", ["C3,4", "C2,4", "C1,4", "C1,5"]),
        
        # Move box4 to goal (C3,1)
        ("box4", ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]),
        
        # Move box5 to goal (C1,3)
        ("box5", ["C4,3", "C3,3", "C2,3", "C1,3"])
    ]
    
    # Execute each move in sequence
    for box, path in moves:
        for new_pos in path:
            if is_valid_move(current_state, box, new_pos):
                current_state = move_box(current_state, box, new_pos)
                # Only add state if it's different from the previous one
                if current_state != solution[-1]:
                    solution.append(current_state.copy())
    
    return solution

# Generate and print solution
solution = generate_solution()
print(json.dumps(solution))