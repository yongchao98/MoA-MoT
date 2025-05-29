import json

initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

def solve():
    solution = []
    current_state = initial_state.copy()
    solution.append(current_state.copy())
    
    # First move box1 to C2,8
    moves_box1 = [
        ("box1", "C1,5"),
        ("box1", "C1,6"),
        ("box1", "C1,7"),
        ("box1", "C1,8"),
        ("box1", "C2,8")
    ]
    
    # Then move box2 to C2,2
    moves_box2 = [
        ("box2", "C2,7"),
        ("box2", "C2,6"),
        ("box2", "C2,5"),
        ("box2", "C2,4"),
        ("box2", "C2,3"),
        ("box2", "C2,2")
    ]
    
    # Then move box3 to C1,2
    moves_box3 = [
        ("box3", "C1,2")
    ]
    
    # Finally move box4 to C2,5
    moves_box4 = [
        ("box4", "C2,2"),
        ("box4", "C2,3"),
        ("box4", "C2,4"),
        ("box4", "C2,5")
    ]
    
    # Execute all moves in sequence
    all_moves = moves_box1 + moves_box2 + moves_box3 + moves_box4
    
    for box, new_pos in all_moves:
        # Check if move is valid (destination is not occupied)
        occupied = set(current_state.values())
        if new_pos not in occupied:
            current_state[box] = new_pos
            solution.append(current_state.copy())
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))