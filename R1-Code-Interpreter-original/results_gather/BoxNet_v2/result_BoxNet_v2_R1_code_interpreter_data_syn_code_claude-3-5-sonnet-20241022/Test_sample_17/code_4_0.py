import json

initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}

goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

def create_solution():
    states = []
    current = initial_state.copy()
    states.append(current.copy())
    
    # Move box1 to goal (C3,6)
    moves = [
        ("box1", "C1,6"),
        ("box1", "C2,6"),
        ("box1", "C3,6"),
        ("box6", "C2,6"),  # Move box6 to goal
        ("box5", "C4,6"),  # Start moving box5
        ("box5", "C4,5"),
        ("box5", "C4,4"),
        ("box5", "C4,3"),
        ("box5", "C3,3"),
        ("box5", "C2,3"),  # box5 at goal
        ("box4", "C3,3"),  # Move box4 to goal
        ("box3", "C3,5"),  # Start moving box3
        ("box3", "C2,5"),
        ("box3", "C2,6"),
        ("box3", "C2,7"),  # box3 at goal
        ("box2", "C3,8"),  # Start moving box2
        ("box2", "C4,8"),
        ("box2", "C4,7"),
        ("box2", "C4,6"),
        ("box2", "C4,5"),
        ("box2", "C4,4"),
        ("box2", "C4,3"),
        ("box2", "C4,2"),
        ("box2", "C4,1")   # box2 at goal
    ]
    
    for box, new_pos in moves:
        current = current.copy()
        current[box] = new_pos
        states.append(current.copy())
    
    return states

# Generate and print solution
solution = create_solution()
print(json.dumps(solution))