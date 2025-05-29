import json

def create_move(state, box, new_pos):
    new_state = state.copy()
    new_state[box] = new_pos
    return new_state

def solve_puzzle():
    solution = []
    
    # Start with initial state
    current_state = {
        "box1": "C5,5",
        "box2": "C1,3",
        "box3": "C3,3",
        "box4": "C4,5",
        "box5": "C4,4"
    }
    solution.append(current_state)
    
    # Move box2 to goal (C3,4)
    moves = [
        ("box2", "C1,4"),
        ("box2", "C2,4"),
        ("box2", "C3,4")
    ]
    for box, pos in moves:
        current_state = create_move(current_state, box, pos)
        solution.append(current_state)
    
    # Move box3 to goal (C1,5)
    moves = [
        ("box3", "C2,3"),
        ("box3", "C2,4"),
        ("box3", "C2,5"),
        ("box3", "C1,5")
    ]
    for box, pos in moves:
        current_state = create_move(current_state, box, pos)
        solution.append(current_state)
    
    # Move box5 to goal (C1,3)
    moves = [
        ("box5", "C4,3"),
        ("box5", "C3,3"),
        ("box5", "C2,3"),
        ("box5", "C1,3")
    ]
    for box, pos in moves:
        current_state = create_move(current_state, box, pos)
        solution.append(current_state)
    
    # Move box4 to goal (C3,1)
    moves = [
        ("box4", "C4,4"),
        ("box4", "C4,3"),
        ("box4", "C4,2"),
        ("box4", "C4,1"),
        ("box4", "C3,1")
    ]
    for box, pos in moves:
        current_state = create_move(current_state, box, pos)
        solution.append(current_state)
    
    # Move box1 to goal (C2,5)
    moves = [
        ("box1", "C4,5"),
        ("box1", "C3,5"),
        ("box1", "C2,5")
    ]
    for box, pos in moves:
        current_state = create_move(current_state, box, pos)
        solution.append(current_state)
    
    return solution

# Generate and print solution
solution = solve_puzzle()
print(json.dumps(solution))