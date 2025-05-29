import json
from collections import deque

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

def get_next_state(current_state, box_to_move):
    occupied = set(current_state.values())
    current_pos = current_state[box_to_move]
    best_next_pos = None
    min_distance = float('inf')
    
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied:
            dist = manhattan_distance(next_pos, goal_state[box_to_move])
            if dist < min_distance:
                min_distance = dist
                best_next_pos = next_pos
    
    if best_next_pos:
        new_state = current_state.copy()
        new_state[box_to_move] = best_next_pos
        return new_state
    return None

def solve():
    solution = [initial_state]
    current_state = initial_state.copy()
    
    # Try to move each box to its goal one by one
    for box in ["box1", "box2", "box3", "box4", "box5", "box6"]:
        while current_state[box] != goal_state[box]:
            next_state = get_next_state(current_state, box)
            if next_state is None:
                # If we can't move this box, try moving other boxes out of the way
                for other_box in current_state:
                    if other_box != box and current_state[other_box] != goal_state[other_box]:
                        temp_state = get_next_state(current_state, other_box)
                        if temp_state:
                            current_state = temp_state
                            solution.append(current_state)
                            break
            else:
                current_state = next_state
                solution.append(current_state)
            
            # Prevent infinite loops
            if len(solution) > 100:
                return None
    
    return solution

# Find and print solution
solution = solve()
if solution:
    # Filter out consecutive duplicate states
    filtered_solution = [solution[0]]
    for state in solution[1:]:
        if state != filtered_solution[-1]:
            filtered_solution.append(state)
    print(json.dumps(filtered_solution))
else:
    print("No solution found")