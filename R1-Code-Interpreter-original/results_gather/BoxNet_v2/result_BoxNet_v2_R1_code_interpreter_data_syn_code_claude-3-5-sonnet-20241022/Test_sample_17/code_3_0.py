import json
from collections import deque

initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}

goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"], "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"], "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"], "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C4,7": ["C4,6", "C4,8", "C3,7"], "C4,8": ["C4,7", "C3,8"]
}

def get_next_moves(state):
    moves = []
    occupied = set(state.values())
    
    # Try to move boxes that are furthest from their goals first
    boxes_to_move = []
    for box in state:
        if state[box] != goal_state[box]:
            boxes_to_move.append(box)
    
    for box in boxes_to_move:
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def solve():
    # Direct solution construction
    solution = [initial_state.copy()]
    current_state = initial_state.copy()
    
    # Move box1 to C3,6
    moves_box1 = ["C1,7", "C1,6", "C2,6", "C3,6"]
    for i in range(len(moves_box1)-1):
        new_state = current_state.copy()
        new_state["box1"] = moves_box1[i+1]
        solution.append(new_state)
        current_state = new_state
    
    # Move box6 to C2,6
    moves_box6 = ["C3,6", "C2,6"]
    for i in range(len(moves_box6)-1):
        new_state = current_state.copy()
        new_state["box6"] = moves_box6[i+1]
        solution.append(new_state)
        current_state = new_state
    
    # Move box5 to C2,3
    moves_box5 = ["C4,7", "C4,6", "C4,5", "C4,4", "C4,3", "C3,3", "C2,3"]
    for i in range(len(moves_box5)-1):
        new_state = current_state.copy()
        new_state["box5"] = moves_box5[i+1]
        solution.append(new_state)
        current_state = new_state
    
    # Move box4 to C3,3
    moves_box4 = ["C2,4", "C3,4", "C3,3"]
    for i in range(len(moves_box4)-1):
        new_state = current_state.copy()
        new_state["box4"] = moves_box4[i+1]
        solution.append(new_state)
        current_state = new_state
    
    # Move box3 to C2,7
    moves_box3 = ["C4,5", "C3,5", "C2,5", "C2,6", "C2,7"]
    for i in range(len(moves_box3)-1):
        new_state = current_state.copy()
        new_state["box3"] = moves_box3[i+1]
        solution.append(new_state)
        current_state = new_state
    
    # Move box2 to C4,1
    moves_box2 = ["C2,8", "C3,8", "C4,8", "C4,7", "C4,6", "C4,5", "C4,4", "C4,3", "C4,2", "C4,1"]
    for i in range(len(moves_box2)-1):
        new_state = current_state.copy()
        new_state["box2"] = moves_box2[i+1]
        solution.append(new_state)
        current_state = new_state
    
    return solution

# Solve and print result
solution = solve()
print(json.dumps(solution))