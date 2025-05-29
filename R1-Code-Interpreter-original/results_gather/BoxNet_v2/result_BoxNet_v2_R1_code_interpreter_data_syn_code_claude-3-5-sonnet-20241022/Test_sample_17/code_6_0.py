import json

def create_solution():
    # Initial state
    current = {
        "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
        "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
    }
    states = [current.copy()]
    
    # Sequence of moves (box, position)
    moves = [
        # First get box1 and box6 to goals
        ("box1", "C2,7"),
        ("box1", "C3,7"),
        ("box1", "C3,6"),
        ("box6", "C2,6"),
        
        # Move box3 to goal
        ("box3", "C3,5"),
        ("box3", "C2,5"),
        ("box3", "C2,6"),
        ("box3", "C2,7"),
        
        # Move box5 to goal
        ("box5", "C4,6"),
        ("box5", "C3,6"),
        ("box5", "C3,5"),
        ("box5", "C3,4"),
        ("box5", "C3,3"),
        ("box5", "C2,3"),
        
        # Move box4 to goal
        ("box4", "C3,4"),
        ("box4", "C3,3"),
        
        # Finally move box2 to goal
        ("box2", "C3,8"),
        ("box2", "C4,8"),
        ("box2", "C4,7"),
        ("box2", "C4,6"),
        ("box2", "C4,5"),
        ("box2", "C4,4"),
        ("box2", "C4,3"),
        ("box2", "C4,2"),
        ("box2", "C4,1")
    ]
    
    # Generate all states
    for box, pos in moves:
        current = current.copy()
        current[box] = pos
        states.append(current)
    
    return states[:10]  # Return only first 10 states to ensure output fits

solution = create_solution()
print(json.dumps(solution))