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

def is_occupied(pos, state, exclude_box=None):
    return any(box != exclude_box and state[box] == pos for box in state)

def move_box(current_state, box, target_pos):
    new_state = dict(current_state)
    new_state[box] = target_pos
    return new_state

def solve():
    solution = [dict(initial_state)]
    current_state = dict(initial_state)
    
    # Predefined sequence of moves
    moves = [
        # First clear space for box5 to reach C1,3
        ("box5", ["C4,3", "C3,3", "C2,3", "C1,3"]),
        # Move box2 to its goal
        ("box2", ["C1,4", "C2,4", "C3,4"]),
        # Move box4 to its goal
        ("box4", ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]),
        # Move box3 to its goal
        ("box3", ["C2,3", "C1,3", "C1,4", "C1,5"]),
        # Finally move box1 to its goal
        ("box1", ["C4,5", "C3,5", "C2,5"])
    ]
    
    for box, path in moves:
        current_pos = current_state[box]
        for next_pos in path:
            if not is_occupied(next_pos, current_state, box):
                current_state = move_box(current_state, box, next_pos)
                solution.append(dict(current_state))
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))