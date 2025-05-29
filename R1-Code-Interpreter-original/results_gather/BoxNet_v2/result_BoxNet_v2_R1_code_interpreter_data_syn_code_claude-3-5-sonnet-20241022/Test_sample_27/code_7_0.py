import json

def make_move(states, box, moves):
    current_state = dict(states[-1])  # Copy the last state
    for move in moves:
        current_state = dict(current_state)  # Make a new copy for each move
        current_state[box] = move
        states.append(current_state)
    return states

def solve():
    # Initialize states list with initial state
    states = [{
        "box1": "C5,5",
        "box2": "C1,3",
        "box3": "C3,3",
        "box4": "C4,5",
        "box5": "C4,4"
    }]
    
    # Define all moves for each box
    moves = [
        # Move box2 to C3,4
        ("box2", ["C1,4", "C2,4", "C3,4"]),
        
        # Move box3 to C1,5
        ("box3", ["C2,3", "C2,4", "C2,5", "C1,5"]),
        
        # Move box5 to C1,3
        ("box5", ["C4,3", "C3,3", "C2,3", "C1,3"]),
        
        # Move box4 to C3,1
        ("box4", ["C4,4", "C4,3", "C4,2", "C4,1", "C3,1"]),
        
        # Move box1 to C2,5
        ("box1", ["C4,5", "C3,5", "C2,5"])
    ]
    
    # Execute all moves sequentially
    for box, box_moves in moves:
        states = make_move(states, box, box_moves)
    
    return states

# Generate solution and convert to JSON string with small chunk size
solution = solve()
json_str = json.dumps(solution)

# Print the solution in smaller chunks to avoid truncation
print(json_str)