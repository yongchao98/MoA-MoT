import json
import heapq
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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
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

def solve_boxes():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
        
        for next_state in get_valid_moves(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
                
                # Early termination if path gets too long
                if len(new_path) > 50:
                    continue
    
    return None

# Solve and format the solution
solution = solve_boxes()
if solution:
    # Ensure the solution is properly formatted
    formatted_solution = []
    for state in solution:
        formatted_solution.append(dict(sorted(state.items())))
    print(json.dumps(formatted_solution))
else:
    print("No solution found")