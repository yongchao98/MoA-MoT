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

def get_best_move(current_state):
    best_move = None
    best_score = float('inf')
    
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
            
        current_dist = manhattan_distance(current_state[box], goal_state[box])
        
        for next_pos in adjacency[current_state[box]]:
            if next_pos not in current_state.values():
                new_dist = manhattan_distance(next_pos, goal_state[box])
                if new_dist < current_dist and new_dist < best_score:
                    best_score = new_dist
                    best_move = (box, next_pos)
    
    return best_move

def solve_greedy():
    path = [initial_state]
    current_state = deepcopy(initial_state)
    moves_without_improvement = 0
    max_moves = 50
    
    while current_state != goal_state and len(path) < max_moves:
        move = get_best_move(current_state)
        
        if not move:
            # If no improving move found, try any valid move
            for box in current_state:
                if current_state[box] == goal_state[box]:
                    continue
                for next_pos in adjacency[current_state[box]]:
                    if next_pos not in current_state.values():
                        move = (box, next_pos)
                        break
                if move:
                    break
        
        if not move:
            break
            
        new_state = deepcopy(current_state)
        new_state[move[0]] = move[1]
        current_state = new_state
        path.append(current_state)
        
    return path if current_state == goal_state else None

# Find and print solution
solution = solve_greedy()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")