import json
from copy import deepcopy

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Simplified adjacency for critical paths
adjacency = {
    "C1,7": ["C2,7", "C1,6"], "C2,7": ["C1,7", "C3,7", "C2,6", "C2,8"],
    "C2,8": ["C2,7", "C3,8"], "C3,8": ["C2,8", "C4,8", "C3,7"],
    "C4,8": ["C3,8", "C4,7"], "C4,7": ["C4,8", "C4,6", "C3,7"],
    "C4,6": ["C4,7", "C4,5", "C3,6"], "C4,5": ["C4,6", "C4,4", "C3,5"],
    "C4,4": ["C4,5", "C4,3", "C3,4"], "C4,3": ["C4,4", "C4,2", "C3,3"],
    "C4,2": ["C4,3", "C4,1", "C3,2"], "C4,1": ["C4,2", "C3,1"],
    "C3,6": ["C3,7", "C3,5", "C2,6", "C4,6"], "C3,5": ["C3,6", "C3,4", "C2,5", "C4,5"],
    "C3,4": ["C3,5", "C3,3", "C2,4", "C4,4"], "C3,3": ["C3,4", "C3,2", "C2,3", "C4,3"],
    "C2,6": ["C2,7", "C2,5", "C1,6", "C3,6"], "C2,5": ["C2,6", "C2,4", "C1,5", "C3,5"],
    "C2,4": ["C2,5", "C2,3", "C1,4", "C3,4"], "C2,3": ["C2,4", "C2,2", "C1,3", "C3,3"]
}

def get_next_moves(state):
    moves = []
    occupied = set(state.values())
    
    # Try to move boxes that are furthest from their goals first
    priorities = []
    for box, pos in state.items():
        if pos != goal_state[box]:
            priorities.append((box, pos))
    
    for box, current_pos in priorities:
        if current_pos in adjacency:
            for next_pos in adjacency[current_pos]:
                if next_pos not in occupied:
                    new_state = deepcopy(state)
                    new_state[box] = next_pos
                    moves.append(new_state)
    
    return moves

def is_closer_to_goal(new_pos, old_pos, goal_pos):
    new_row, new_col = int(new_pos.split(',')[0][1]), int(new_pos.split(',')[1])
    old_row, old_col = int(old_pos.split(',')[0][1]), int(old_pos.split(',')[1])
    goal_row, goal_col = int(goal_pos.split(',')[0][1]), int(goal_pos.split(',')[1])
    
    old_dist = abs(old_row - goal_row) + abs(old_col - goal_col)
    new_dist = abs(new_row - goal_row) + abs(new_col - goal_col)
    return new_dist < old_dist

def solve():
    current_state = initial_state
    solution = [current_state]
    visited = {json.dumps(current_state, sort_keys=True)}
    
    while current_state != goal_state:
        best_next_state = None
        best_improvement = -1
        
        for next_state in get_next_moves(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str in visited:
                continue
                
            improvement = 0
            for box in next_state:
                if next_state[box] != current_state[box]:
                    if is_closer_to_goal(next_state[box], current_state[box], goal_state[box]):
                        improvement += 1
            
            if improvement > best_improvement:
                best_improvement = improvement
                best_next_state = next_state
        
        if best_next_state is None:
            return None
            
        current_state = best_next_state
        solution.append(current_state)
        visited.add(json.dumps(current_state, sort_keys=True))
        
        if len(solution) > 30:  # Limit solution length
            return None
    
    return solution

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")