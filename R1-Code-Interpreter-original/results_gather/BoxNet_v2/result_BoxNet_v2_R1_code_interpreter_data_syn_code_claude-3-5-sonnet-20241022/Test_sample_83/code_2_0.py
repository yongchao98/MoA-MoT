import json
from copy import deepcopy
from heapq import heappush, heappop
from collections import deque

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

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box, pos in state.items():
        if pos == goal_state[box]:
            continue
            
        for next_pos in adjacency[pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def total_distance(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def solve():
    start = initial_state
    visited = {state_to_string(start)}
    queue = deque([(start, [start])])
    best_distance = float('inf')
    best_path = None
    
    while queue:
        current, path = queue.popleft()
        current_distance = total_distance(current)
        
        if current_distance < best_distance:
            best_distance = current_distance
            best_path = path
            
        if current == goal_state:
            return path
            
        for next_state in get_valid_moves(current):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
                
        if len(queue) > 10000:  # Prevent memory overflow
            queue = deque(sorted(queue, key=lambda x: total_distance(x[0]))[:5000])
    
    return best_path

# Solve and print the result
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")