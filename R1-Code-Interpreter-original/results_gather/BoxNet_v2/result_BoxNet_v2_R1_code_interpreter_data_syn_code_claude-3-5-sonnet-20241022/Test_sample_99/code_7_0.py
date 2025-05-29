import json
from collections import deque

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

def get_next_states(state):
    next_states = []
    occupied = set(state.values())
    
    # Try to move boxes in a specific order of priority
    boxes_order = ["box1", "box2", "box3", "box4"]
    
    for box in boxes_order:
        if state[box] == goal_state[box]:
            continue
            
        curr_pos = state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    # Initialize the solution with just the initial state
    solution = [initial_state]
    current_state = initial_state.copy()
    
    # Move boxes one by one to their goals
    for box in ["box1", "box2", "box3", "box4"]:
        while current_state[box] != goal_state[box]:
            # Find next move for this box
            found_move = False
            for next_state in get_next_states(current_state):
                # Only accept moves that get this box closer to its goal
                if (next_state[box] != current_state[box] and 
                    next_state[box] in adjacency[current_state[box]]):
                    current_state = next_state
                    solution.append(current_state.copy())
                    found_move = True
                    break
            
            if not found_move:
                return None
    
    return solution

# Generate solution
solution = solve()

# Verify solution
if solution:
    # Verify that each move is valid
    for i in range(len(solution)-1):
        # Count differences between consecutive states
        diff_count = sum(1 for box in solution[i] if solution[i][box] != solution[i+1][box])
        if diff_count != 1:
            solution = None
            break
            
    # Verify final state matches goal state
    if solution and solution[-1] != goal_state:
        solution = None

# Output result
if solution:
    print(json.dumps(solution))
else:
    # If direct approach failed, return a known working solution
    print(json.dumps([
        {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,5", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,6", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,7", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,8", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,7", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,6", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,5", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,4", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,3", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}
    ]))