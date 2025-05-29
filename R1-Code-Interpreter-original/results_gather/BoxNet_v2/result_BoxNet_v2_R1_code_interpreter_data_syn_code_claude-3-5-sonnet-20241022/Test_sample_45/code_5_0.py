import json
from copy import deepcopy

initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Simplified adjacency for faster lookup
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

def solve():
    solution = [initial_state]
    current_state = dict(initial_state)
    
    # Move boxes one by one to their goals
    move_sequence = [
        ("box1", "C3,6"),  # Move box1 to its goal
        ("box6", "C2,6"),  # Move box6 to its goal
        ("box5", "C2,3"),  # Move box5 to its goal
        ("box4", "C3,3"),  # Move box4 to its goal
        ("box3", "C2,7"),  # Move box3 to its goal
        ("box2", "C4,1")   # Move box2 to its goal
    ]
    
    for box, target in move_sequence:
        while current_state[box] != target:
            # Find next move towards target
            current_pos = current_state[box]
            occupied = set(current_state.values())
            
            # Get best next position
            best_next = None
            min_dist = float('inf')
            
            for next_pos in adjacency[current_pos]:
                if next_pos not in occupied:
                    # Simple distance calculation
                    curr_row, curr_col = int(next_pos.split(',')[0][1]), int(next_pos.split(',')[1])
                    target_row, target_col = int(target.split(',')[0][1]), int(target.split(',')[1])
                    dist = abs(curr_row - target_row) + abs(curr_col - target_col)
                    
                    if dist < min_dist:
                        min_dist = dist
                        best_next = next_pos
            
            if best_next:
                new_state = deepcopy(current_state)
                new_state[box] = best_next
                current_state = new_state
                solution.append(new_state)
            else:
                # If stuck, move other boxes out of the way
                for other_box in current_state:
                    if other_box != box and current_state[other_box] != goal_state[other_box]:
                        for next_pos in adjacency[current_state[other_box]]:
                            if next_pos not in occupied:
                                new_state = deepcopy(current_state)
                                new_state[other_box] = next_pos
                                current_state = new_state
                                solution.append(new_state)
                                break
                        break
    
    return solution

# Find and print solution
solution = solve()
print(json.dumps(solution))