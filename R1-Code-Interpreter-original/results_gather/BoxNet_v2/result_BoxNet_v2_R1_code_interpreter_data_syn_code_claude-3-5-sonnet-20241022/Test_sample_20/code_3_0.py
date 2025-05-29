import json
from collections import deque
from copy import deepcopy

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_valid_moves(state, adjacency, goal_state):
    moves = []
    occupied = set(state.values())
    
    # Try to move boxes that are closest to their goals first
    for box, current_pos in state.items():
        if current_pos == goal_state[box]:
            continue
            
        # Calculate distance to goal
        dist_to_goal = manhattan_distance(current_pos, goal_state[box])
        
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                # Calculate if this move gets us closer to the goal
                new_dist = manhattan_distance(next_pos, goal_state[box])
                if new_dist < dist_to_goal:  # Only consider moves that get closer to goal
                    new_state = deepcopy(state)
                    new_state[box] = next_pos
                    moves.append(new_state)
                    
                    # If this move reaches the goal, prioritize it
                    if next_pos == goal_state[box]:
                        return [new_state]
    
    return moves

def solve_puzzle():
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
        "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
        "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
        "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
        "C2,5": ["C2,4", "C1,5", "C3,5"], "C3,1": ["C3,2", "C2,1", "C4,1"],
        "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
        "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
        "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
        "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
        "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
        "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
        "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
        "C4,5": ["C4,4", "C3,5", "C5,5"], "C5,1": ["C5,2", "C4,1"],
        "C5,2": ["C5,1", "C5,3", "C4,2"], "C5,3": ["C5,2", "C5,4", "C4,3"],
        "C5,4": ["C5,3", "C5,5", "C4,4"], "C5,5": ["C5,4", "C4,5"]
    }

    queue = deque([(initial_state, [initial_state])])
    visited = {json.dumps(initial_state, sort_keys=True)}
    max_depth = 20  # Limit search depth
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
            
        if len(path) > max_depth:
            continue
            
        for next_state in get_valid_moves(current_state, adjacency, goal_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
    
    return None

# Find and print solution
solution = solve_puzzle()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")