import json

def get_next_states(state):
    moves = []
    occupied = set(state.values())
    
    # Only move one box that's not in its goal position
    for box in state:
        if state[box] == goal_state[box]:
            continue
            
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def direct_solve():
    # Initial state
    path = [{"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"}]
    
    # First move box1 to C2,8
    moves_box1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for pos in moves_box1:
        new_state = path[-1].copy()
        new_state["box1"] = pos
        path.append(new_state)
    
    # Then move box2 to C2,2
    moves_box2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for pos in moves_box2:
        new_state = path[-1].copy()
        new_state["box2"] = pos
        path.append(new_state)
    
    # Then move box3 to C1,2
    new_state = path[-1].copy()
    new_state["box3"] = "C1,2"
    path.append(new_state)
    
    # Finally move box4 to C2,5
    moves_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for pos in moves_box4:
        new_state = path[-1].copy()
        new_state["box4"] = pos
        path.append(new_state)
    
    return path

# Get solution
solution = direct_solve()
print(json.dumps(solution))