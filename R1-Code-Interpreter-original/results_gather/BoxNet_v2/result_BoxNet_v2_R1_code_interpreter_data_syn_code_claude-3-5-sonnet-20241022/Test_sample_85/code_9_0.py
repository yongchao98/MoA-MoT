import json
from collections import deque
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

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def solve():
    # Start with moving box4 out of the way
    solution = [
        {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},  # Initial state
        {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,3"},  # Move box4 to C2,3
        {"box1": "C1,4", "box2": "C2,8", "box3": "C1,2", "box4": "C2,3"},  # Move box3 to C1,2 (goal)
        {"box1": "C1,4", "box2": "C2,7", "box3": "C1,2", "box4": "C2,3"},  # Move box2 towards goal
        {"box1": "C1,4", "box2": "C2,6", "box3": "C1,2", "box4": "C2,3"},  # Move box2
        {"box1": "C1,4", "box2": "C2,5", "box3": "C1,2", "box4": "C2,3"},  # Move box2
        {"box1": "C1,4", "box2": "C2,4", "box3": "C1,2", "box4": "C2,3"},  # Move box2
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},  # Move box2 to goal
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},  # Move box4 towards goal
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # Move box4 to goal
        {"box1": "C1,5", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # Move box1 towards goal
        {"box1": "C1,6", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # Move box1
        {"box1": "C1,7", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # Move box1
        {"box1": "C1,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},  # Move box1
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}   # Move box1 to goal
    ]
    
    # Verify each move is valid
    for i in range(len(solution) - 1):
        current_state = solution[i]
        next_state = solution[i + 1]
        
        # Check if the move is valid
        valid_moves = get_valid_moves(current_state)
        if next_state not in valid_moves:
            return None
    
    # Verify final state matches goal state
    if solution[-1] != goal_state:
        return None
        
    return solution

# Solve and print result
solution = solve()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")