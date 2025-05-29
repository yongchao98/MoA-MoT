import heapq
from copy import deepcopy
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

def get_valid_moves(state):
    moves = []
    for box in state:
        if state[box] == goal_state[box]:
            continue
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def solve():
    queue = [(initial_state, [initial_state])]
    visited = {json.dumps(initial_state, sort_keys=True)}
    
    while queue:
        current_state, path = queue.pop(0)
        
        if current_state == goal_state:
            return path
            
        for next_state in get_valid_moves(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

# Solve and print result
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")