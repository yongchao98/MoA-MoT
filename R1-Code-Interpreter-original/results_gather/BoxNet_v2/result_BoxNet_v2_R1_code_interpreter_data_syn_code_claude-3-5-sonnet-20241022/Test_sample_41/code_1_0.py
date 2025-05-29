import json
import heapq
from collections import defaultdict

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

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        for next_pos in adjacency[current_state[box]]:
            # Check if next_pos is occupied by another box
            if next_pos in current_state.values():
                continue
            new_state = current_state.copy()
            new_state[box] = next_pos
            next_states.append(new_state)
    return next_states

def solve():
    start = initial_state
    visited = {state_to_string(start)}
    pq = [(heuristic(start), 0, start, [start])]
    max_depth = 30
    
    while pq:
        _, depth, current, path = heapq.heappop(pq)
        
        if current == goal_state:
            return path
            
        if depth >= max_depth:
            continue
            
        for next_state in get_next_states(current):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h_score = heuristic(next_state)
                heapq.heappush(pq, (h_score + depth + 1, depth + 1, next_state, path + [next_state]))
    
    return None

solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")