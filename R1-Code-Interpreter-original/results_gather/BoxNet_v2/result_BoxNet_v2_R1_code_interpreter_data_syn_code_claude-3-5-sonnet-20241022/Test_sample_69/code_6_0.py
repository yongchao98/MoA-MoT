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

def is_position_free(pos, state, moving_box):
    return not any(box != moving_box and state[box] == pos for box in state)

def generate_solution():
    solution = [dict(initial_state)]
    current = dict(initial_state)
    
    # Move box1 to goal
    moves = [("box1", ["C4,5", "C3,5", "C2,5"])]
    
    # Move box5 out of the way first
    moves.extend([
        ("box5", ["C3,4", "C2,4", "C1,4"]),
        ("box2", ["C2,3", "C3,3", "C3,4"]),
        ("box4", ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]),
        ("box3", ["C2,3", "C1,3", "C1,4", "C1,5"]),
        ("box5", ["C1,3"])
    ])
    
    for box, path in moves:
        for next_pos in path:
            if is_position_free(next_pos, current, box):
                current = move_box(current, box, next_pos)
                if current != solution[-1]:  # Only add if state has changed
                    solution.append(dict(current))
    
    return solution

# Generate and print solution
solution = generate_solution()
print(json.dumps(solution))