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

def is_valid_move(current_pos, next_pos, state, box):
    if next_pos not in adjacency[current_pos]:
        return False
    return not any(b != box and state[b] == next_pos for b in state)

def move_box(state, box, path):
    result = []
    current = dict(state)
    
    for next_pos in path:
        if is_valid_move(current[box], next_pos, current, box):
            current = dict(current)
            current[box] = next_pos
            result.append(current)
    
    return result

def solve():
    solution = [initial_state]
    current = dict(initial_state)
    
    # Carefully ordered sequence of moves
    moves = [
        ("box1", ["C3,5", "C2,5"]),  # Move box1 to goal
        ("box5", ["C3,4", "C2,4", "C1,4", "C1,3"]),  # Move box5 to goal
        ("box2", ["C2,3", "C3,3", "C3,4"]),  # Move box2 to goal
        ("box4", ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]),  # Move box4 to goal
        ("box3", ["C2,3", "C1,3", "C1,4", "C1,5"])  # Move box3 to goal
    ]
    
    for box, path in moves:
        if current[box] != path[-1]:  # Only move if not already at destination
            new_states = move_box(current, box, path)
            for state in new_states:
                current = state
                solution.append(current)
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))