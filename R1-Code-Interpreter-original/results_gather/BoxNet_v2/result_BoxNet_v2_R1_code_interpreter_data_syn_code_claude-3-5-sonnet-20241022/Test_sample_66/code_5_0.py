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

def get_box_moves(state, box):
    moves = []
    occupied = set(state.values())
    current_pos = state[box]
    target_pos = goal_state[box]
    
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied:
            new_state = deepcopy(state)
            new_state[box] = next_pos
            # Calculate if this move gets us closer to the goal
            old_dist = manhattan_distance(current_pos, target_pos)
            new_dist = manhattan_distance(next_pos, target_pos)
            moves.append((new_dist - old_dist, new_state))
    
    moves.sort()  # Sort by distance improvement
    return [m[1] for m in moves]

def solve_for_box(state, box, depth=0, max_depth=5, visited=None):
    if visited is None:
        visited = set()
    
    if depth > max_depth:
        return None
    
    if state[box] == goal_state[box]:
        return [state]
    
    state_str = json.dumps(state, sort_keys=True)
    if state_str in visited:
        return None
    visited.add(state_str)
    
    for next_state in get_box_moves(state, box):
        result = solve_for_box(next_state, box, depth + 1, max_depth, visited)
        if result is not None:
            return [state] + result
    
    return None

def solve():
    current_state = initial_state
    path = [current_state]
    boxes_to_move = list(initial_state.keys())
    visited_states = set()
    
    # Try to solve for each box one at a time
    while boxes_to_move:
        # Find the box that's furthest from its goal
        box_distances = [(manhattan_distance(current_state[box], goal_state[box]), box) 
                        for box in boxes_to_move]
        box_distances.sort(reverse=True)
        
        moved = False
        for _, box in box_distances:
            if current_state[box] == goal_state[box]:
                boxes_to_move.remove(box)
                moved = True
                continue
                
            # Try to move this box closer to its goal
            result = solve_for_box(current_state, box)
            if result is not None:
                # Update current state and path
                current_state = result[-1]
                path.extend(result[1:])
                moved = True
                if current_state[box] == goal_state[box]:
                    boxes_to_move.remove(box)
                break
        
        if not moved:
            return None
        
        # Check if we're in a cycle
        state_str = json.dumps(current_state, sort_keys=True)
        if state_str in visited_states:
            return None
        visited_states.add(state_str)
        
        # Check if we've reached the goal
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
    
    return path if all(current_state[box] == goal_state[box] for box in current_state) else None

# Find and print the solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")