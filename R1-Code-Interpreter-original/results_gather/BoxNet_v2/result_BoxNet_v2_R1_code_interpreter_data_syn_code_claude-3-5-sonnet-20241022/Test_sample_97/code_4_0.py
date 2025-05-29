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

def get_priority_box(state):
    distances = []
    for box in state:
        if state[box] != goal_state[box]:
            dist = manhattan_distance(state[box], goal_state[box])
            distances.append((dist, box))
    return sorted(distances, reverse=True)[0][1] if distances else None

def get_next_moves(state, box):
    moves = []
    current_pos = state[box]
    for next_pos in adjacency[current_pos]:
        if next_pos not in state.values():
            new_state = deepcopy(state)
            new_state[box] = next_pos
            moves.append(new_state)
    return moves

def solve_step(current_state, path, depth, max_depth, visited):
    if depth > max_depth:
        return None
    
    if current_state == goal_state:
        return path
    
    priority_box = get_priority_box(current_state)
    if not priority_box:
        return path
    
    next_states = get_next_moves(current_state, priority_box)
    for next_state in next_states:
        state_str = str(sorted(next_state.items()))
        if state_str not in visited:
            visited.add(state_str)
            result = solve_step(next_state, path + [next_state], depth + 1, max_depth, visited)
            if result:
                return result
    
    return None

def solve():
    max_depth = 15  # Limit the search depth
    visited = {str(sorted(initial_state.items()))}
    solution = solve_step(initial_state, [initial_state], 0, max_depth, visited)
    return solution

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")