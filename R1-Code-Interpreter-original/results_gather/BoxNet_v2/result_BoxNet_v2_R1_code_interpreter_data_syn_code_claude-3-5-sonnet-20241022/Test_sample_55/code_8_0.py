import json
from copy import deepcopy
import heapq

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

def get_single_move(state, box):
    current_pos = state[box]
    target_pos = goal_state[box]
    
    if current_pos == target_pos:
        return None
        
    occupied = set(state.values())
    best_move = None
    best_dist = float('inf')
    
    for next_pos in adjacency[current_pos]:
        if next_pos not in occupied:
            dist = manhattan_distance(next_pos, target_pos)
            if dist < best_dist:
                best_dist = dist
                best_move = next_pos
                
    return best_move

def solve_greedy():
    current_state = deepcopy(initial_state)
    solution = [current_state]
    moves_count = 0
    max_moves = 100
    
    while current_state != goal_state and moves_count < max_moves:
        # Find box that's furthest from its goal
        boxes = list(current_state.keys())
        boxes.sort(key=lambda b: -manhattan_distance(current_state[b], goal_state[b]))
        
        moved = False
        for box in boxes:
            if current_state[box] != goal_state[box]:
                next_pos = get_single_move(current_state, box)
                if next_pos:
                    new_state = deepcopy(current_state)
                    new_state[box] = next_pos
                    current_state = new_state
                    solution.append(current_state)
                    moved = True
                    moves_count += 1
                    break
        
        if not moved:
            # Try moving any box that's blocking the path
            for box in boxes:
                if current_state[box] != goal_state[box]:
                    for next_pos in adjacency[current_state[box]]:
                        if next_pos not in current_state.values():
                            new_state = deepcopy(current_state)
                            new_state[box] = next_pos
                            current_state = new_state
                            solution.append(current_state)
                            moved = True
                            moves_count += 1
                            break
                if moved:
                    break
        
        if not moved:
            return None
    
    return solution if current_state == goal_state else None

# Find solution
solution = solve_greedy()

# Verify and format solution
if solution:
    # Verify moves are valid
    valid = True
    for i in range(len(solution) - 1):
        diff_count = sum(1 for k in solution[i] if solution[i][k] != solution[i+1][k])
        if diff_count != 1:
            valid = False
            break
    
    if valid:
        formatted_solution = []
        for state in solution:
            formatted_solution.append(dict(sorted(state.items())))
        print(json.dumps(formatted_solution))
    else:
        print("No valid solution found")
else:
    print("No valid solution found")