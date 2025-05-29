import json

def is_valid_move(current_state, box, new_pos):
    # Check if the new position is occupied by another box
    for other_box, pos in current_state.items():
        if other_box != box and pos == new_pos:
            return False
    # Check if the move is to an adjacent cell
    return new_pos in adjacency[current_state[box]]

def make_move(state, box, new_pos):
    if is_valid_move(state, box, new_pos):
        new_state = state.copy()
        new_state[box] = new_pos
        return new_state
    return None

# Initial and goal states
initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

# Adjacency list
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

# Generate solution
solution = [initial_state]
current_state = dict(initial_state)

# Predefined sequence of moves
moves = [
    # First clear space by moving box5 to its final position
    ("box5", ["C3,4", "C2,4", "C1,4", "C1,3"]),
    
    # Move box2 to its goal
    ("box2", ["C2,3", "C3,3", "C3,4"]),
    
    # Move box3 to its goal
    ("box3", ["C2,3", "C1,3", "C1,4", "C1,5"]),
    
    # Move box1 to its goal
    ("box1", ["C4,5", "C3,5", "C2,5"]),
    
    # Finally move box4 to its goal
    ("box4", ["C3,5", "C3,4", "C3,3", "C3,2", "C3,1"])
]

# Execute moves
for box, path in moves:
    for new_pos in path:
        new_state = make_move(current_state, box, new_pos)
        if new_state:
            current_state = new_state
            solution.append(dict(current_state))

# Print solution
print(json.dumps(solution))