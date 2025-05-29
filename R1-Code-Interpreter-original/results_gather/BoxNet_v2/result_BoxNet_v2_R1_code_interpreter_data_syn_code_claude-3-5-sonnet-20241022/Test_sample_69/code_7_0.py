import json

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

def move_box(state, box, new_pos):
    new_state = dict(state)
    new_state[box] = new_pos
    return new_state

def generate_solution():
    solution = []
    current = dict(initial_state)
    solution.append(current)
    
    # Fixed sequence of moves
    moves = [
        # Move box1 to goal
        ("box1", ["C3,5", "C2,5"]),
        # Move box5 to final position
        ("box5", ["C3,4", "C2,4", "C1,4", "C1,3"]),
        # Move box2 to goal
        ("box2", ["C2,3", "C3,3", "C3,4"]),
        # Move box4 to goal
        ("box4", ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]),
        # Move box3 to goal
        ("box3", ["C2,3", "C1,3", "C1,4", "C1,5"])
    ]
    
    for box, path in moves:
        for pos in path:
            current = move_box(current, box, pos)
            solution.append(dict(current))
    
    return solution

# Generate and print solution
solution = generate_solution()
print(json.dumps(solution))