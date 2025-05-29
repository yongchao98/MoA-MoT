import json
import heapq
from copy import deepcopy

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

# Simplified adjacency for clarity
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

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        # Skip if box is at goal
        if state[box] == goal_state[box]:
            continue
            
        curr_pos = state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied or next_pos == state[box]:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def solve_with_iterative_deepening():
    max_depth = 40
    visited = set()
    
    def dfs(current_state, depth, path):
        if depth > max_depth:
            return None
            
        if current_state == goal_state:
            return path
            
        state_str = state_to_string(current_state)
        if state_str in visited:
            return None
            
        visited.add(state_str)
        
        # Sort moves by heuristic value
        moves = get_valid_moves(current_state)
        moves.sort(key=lambda x: sum(manhattan_distance(x[box], goal_state[box]) 
                                   for box in x if x[box] != goal_state[box]))
        
        for next_state in moves:
            # Verify only one box moved
            differences = sum(1 for box in current_state 
                           if current_state[box] != next_state[box])
            if differences != 1:
                continue
                
            result = dfs(next_state, depth + 1, path + [next_state])
            if result:
                return result
                
        visited.remove(state_str)
        return None
    
    return dfs(initial_state, 0, [initial_state])

# Find solution
solution = solve_with_iterative_deepening()

# Format and print solution
if solution:
    # Verify solution
    valid = True
    for i in range(len(solution) - 1):
        differences = sum(1 for box in solution[i] 
                        if solution[i][box] != solution[i+1][box])
        if differences != 1:
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
    print("No solution found")