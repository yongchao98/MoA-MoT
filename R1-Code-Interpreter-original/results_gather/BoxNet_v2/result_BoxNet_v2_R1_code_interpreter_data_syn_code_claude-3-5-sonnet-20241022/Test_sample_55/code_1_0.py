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

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def is_valid_move(current_state, next_state):
    # Check only one box moved
    differences = sum(1 for box in current_state if current_state[box] != next_state[box])
    if differences != 1:
        return False
    
    # Check the move is to an adjacent cell
    for box in current_state:
        if current_state[box] != next_state[box]:
            if next_state[box] not in adjacency[current_state[box]]:
                return False
            # Check destination is not occupied by another box
            for other_box in current_state:
                if other_box != box and next_state[box] == next_state[other_box]:
                    return False
    return True

def generate_solution():
    solution = [initial_state]
    current_state = deepcopy(initial_state)
    
    # Predefined sequence of moves
    moves = [
        # First clear space around target positions
        ("box5", "C4,3"),  # Move box5 away from box4
        ("box5", "C3,3"),  # Continue moving box5
        ("box4", "C4,4"),  # Move box4 to make way
        ("box4", "C3,4"),  # Continue moving box4
        ("box4", "C3,3"),  # Move box4 towards its goal
        ("box4", "C3,2"),  # Continue moving box4
        ("box4", "C3,1"),  # Box4 reaches its goal
        
        # Move box2 to its goal
        ("box2", "C1,4"),
        ("box2", "C2,4"),
        ("box2", "C3,4"),  # Box2 reaches its goal
        
        # Move box5 to its final position
        ("box5", "C2,3"),
        ("box5", "C1,3"),  # Box5 reaches its goal
        
        # Move box3 to its goal
        ("box3", "C2,3"),
        ("box3", "C1,3"),
        ("box3", "C1,4"),
        ("box3", "C1,5"),  # Box3 reaches its goal
        
        # Finally move box1 to its goal
        ("box1", "C4,5"),
        ("box1", "C3,5"),
        ("box1", "C2,5"),  # Box1 reaches its goal
    ]
    
    # Execute moves
    for box, target in moves:
        if current_state[box] == target:
            continue
            
        new_state = deepcopy(current_state)
        new_state[box] = target
        
        if is_valid_move(current_state, new_state):
            current_state = new_state
            solution.append(current_state)
        else:
            return None
    
    # Verify final state matches goal state
    if current_state == goal_state:
        return solution
    return None

# Generate and verify solution
solution = generate_solution()

if solution:
    # Verify all moves are valid
    valid = True
    for i in range(len(solution) - 1):
        if not is_valid_move(solution[i], solution[i+1]):
            valid = False
            break
    
    if valid and solution[-1] == goal_state:
        formatted_solution = []
        for state in solution:
            formatted_solution.append(dict(sorted(state.items())))
        print(json.dumps(formatted_solution))
    else:
        print("No valid solution found")
else:
    print("No valid solution found")