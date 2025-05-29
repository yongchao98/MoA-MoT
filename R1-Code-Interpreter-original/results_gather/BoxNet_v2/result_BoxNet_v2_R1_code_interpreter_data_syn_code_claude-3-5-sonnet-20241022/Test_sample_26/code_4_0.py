import json
from copy import deepcopy

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

def get_valid_moves(state, box):
    moves = []
    current_pos = state[box]
    occupied = set(state.values())
    
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied:
            new_state = deepcopy(state)
            new_state[box] = next_pos
            moves.append(new_state)
    
    return moves

def solve():
    # Create solution path starting with initial state
    solution = [initial_state]
    current_state = initial_state
    
    # Move box1 to C2,8
    while current_state["box1"] != "C2,8":
        moves = get_valid_moves(current_state, "box1")
        for move in moves:
            if move["box1"] in ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]:
                current_state = move
                solution.append(current_state)
                break
    
    # Move box2 to C2,2
    while current_state["box2"] != "C2,2":
        moves = get_valid_moves(current_state, "box2")
        for move in moves:
            if move["box2"] in ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]:
                current_state = move
                solution.append(current_state)
                break
    
    # Move box3 to C1,2
    while current_state["box3"] != "C1,2":
        moves = get_valid_moves(current_state, "box3")
        for move in moves:
            if move["box3"] == "C1,2":
                current_state = move
                solution.append(current_state)
                break
    
    # Move box4 to C2,5
    while current_state["box4"] != "C2,5":
        moves = get_valid_moves(current_state, "box4")
        for move in moves:
            if move["box4"] in ["C2,2", "C2,3", "C2,4", "C2,5"]:
                current_state = move
                solution.append(current_state)
                break
    
    return solution

# Solve and print result
solution = solve()
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")