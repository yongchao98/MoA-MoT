import json
from collections import deque
import copy

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

def is_valid_move(current_state, box, new_pos):
    for other_box, pos in current_state.items():
        if other_box != box and pos == new_pos:
            return False
    return True

def get_manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_total_distance(state):
    return sum(get_manhattan_distance(state[box], goal_state[box]) for box in state)

def solve_for_box(current_state, target_box):
    if current_state[target_box] == goal_state[target_box]:
        return [current_state]
    
    queue = deque([(current_state, [current_state])])
    visited = {str(sorted(current_state.items()))}
    
    while queue:
        state, path = queue.popleft()
        current_pos = state[target_box]
        
        for next_pos in adjacency[current_pos]:
            if is_valid_move(state, target_box, next_pos):
                new_state = copy.deepcopy(state)
                new_state[target_box] = next_pos
                state_str = str(sorted(new_state.items()))
                
                if state_str not in visited:
                    visited.add(state_str)
                    new_path = path + [new_state]
                    
                    if new_state[target_box] == goal_state[target_box]:
                        return new_path
                    
                    queue.append((new_state, new_path))
    return None

def solve():
    # Order boxes by priority (distance to goal and dependencies)
    current_state = copy.deepcopy(initial_state)
    solution = [current_state]
    
    # Move box1 to goal first (it needs to be at C2,8)
    path = solve_for_box(current_state, "box1")
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
    
    # Move box2 to goal (C2,2)
    path = solve_for_box(current_state, "box2")
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
    
    # Move box3 to goal (C1,2)
    path = solve_for_box(current_state, "box3")
    if path:
        solution.extend(path[1:])
        current_state = path[-1]
    
    # Finally move box4 to goal (C2,5)
    path = solve_for_box(current_state, "box4")
    if path:
        solution.extend(path[1:])
    
    return solution

solution = solve()
print(json.dumps(solution))