import json
from copy import deepcopy

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

# Simplified adjacency for direct paths
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

def is_valid_move(state, box, new_pos):
    return new_pos in adjacency[state[box]] and new_pos not in state.values()

def try_move_box(state, box, target):
    if state[box] == target:
        return state, True
    
    current = state[box]
    current_x, current_y = map(int, current.replace('C', '').split(','))
    target_x, target_y = map(int, target.replace('C', '').split(','))
    
    # Try moving horizontally first, then vertically
    new_x = current_x + (1 if target_x > current_x else -1 if target_x < current_x else 0)
    new_y = current_y
    
    new_pos = f"C{new_x},{new_y}"
    if new_x >= 1 and new_x <= 5 and new_pos in adjacency[current] and is_valid_move(state, box, new_pos):
        new_state = deepcopy(state)
        new_state[box] = new_pos
        return new_state, False
    
    # Try moving vertically
    new_x = current_x
    new_y = current_y + (1 if target_y > current_y else -1 if target_y < current_y else 0)
    
    new_pos = f"C{new_x},{new_y}"
    if new_y >= 1 and new_y <= 5 and new_pos in adjacency[current] and is_valid_move(state, box, new_pos):
        new_state = deepcopy(state)
        new_state[box] = new_pos
        return new_state, False
    
    return state, False

def solve():
    current_state = deepcopy(initial_state)
    solution = [current_state]
    
    # Move boxes one by one to their goals
    move_sequence = ["box4", "box5", "box2", "box3", "box1"]
    
    for box in move_sequence:
        while current_state[box] != goal_state[box]:
            new_state, done = try_move_box(current_state, box, goal_state[box])
            if new_state == current_state:
                # Try moving other boxes temporarily to make way
                for other_box in current_state:
                    if other_box != box and current_state[other_box] != goal_state[other_box]:
                        temp_state, _ = try_move_box(current_state, other_box, "C3,3")
                        if temp_state != current_state:
                            current_state = temp_state
                            solution.append(current_state)
                            break
                else:
                    return None
            else:
                current_state = new_state
                solution.append(current_state)
            
            if len(solution) > 50:  # Prevent infinite loops
                return None
    
    return solution

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")