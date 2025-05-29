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

# Adjacency dictionary
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

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_best_move(state, box):
    current_pos = state[box]
    target_pos = goal_state[box]
    current_dist = manhattan_distance(current_pos, target_pos)
    
    best_move = None
    best_dist = current_dist
    occupied = set(state.values())
    
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied:
            dist = manhattan_distance(next_pos, target_pos)
            if dist < best_dist:
                best_dist = dist
                best_move = next_pos
                
    return best_move

def solve():
    current_state = deepcopy(initial_state)
    path = [current_state]
    visited = {json.dumps(current_state, sort_keys=True)}
    max_steps = 100
    
    while len(path) < max_steps:
        # Find box that needs to move
        boxes_to_move = []
        for box in current_state:
            if current_state[box] != goal_state[box]:
                dist = manhattan_distance(current_state[box], goal_state[box])
                boxes_to_move.append((dist, box))
        
        if not boxes_to_move:
            return path
            
        # Sort boxes by distance to goal
        boxes_to_move.sort(reverse=True)
        
        # Try to move each box
        moved = False
        for _, box in boxes_to_move:
            best_move = get_best_move(current_state, box)
            if best_move:
                new_state = deepcopy(current_state)
                new_state[box] = best_move
                state_str = json.dumps(new_state, sort_keys=True)
                
                if state_str not in visited:
                    current_state = new_state
                    path.append(current_state)
                    visited.add(state_str)
                    moved = True
                    break
        
        if not moved:
            # If no box can move closer to its goal, try any valid move
            for _, box in boxes_to_move:
                occupied = set(current_state.values())
                current_pos = current_state[box]
                
                for next_pos in adjacency[current_pos]:
                    if next_pos not in occupied:
                        new_state = deepcopy(current_state)
                        new_state[box] = next_pos
                        state_str = json.dumps(new_state, sort_keys=True)
                        
                        if state_str not in visited:
                            current_state = new_state
                            path.append(current_state)
                            visited.add(state_str)
                            moved = True
                            break
                
                if moved:
                    break
        
        if not moved:
            return None
    
    return None

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")