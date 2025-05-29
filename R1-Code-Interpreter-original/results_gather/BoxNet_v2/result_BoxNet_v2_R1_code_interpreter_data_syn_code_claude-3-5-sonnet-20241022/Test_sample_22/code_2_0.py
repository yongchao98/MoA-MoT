import heapq
from copy import deepcopy
import json

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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def get_next_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
        
        curr_pos = state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append((box, next_pos, new_state))
    
    return moves

def solve_with_bfs():
    from collections import deque
    
    start = json.dumps(initial_state)
    queue = deque([(initial_state, [initial_state])])
    visited = {start}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
            
        for box, next_pos, new_state in get_next_moves(current_state):
            state_str = json.dumps(new_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((new_state, path + [new_state]))
    
    return None

# Solve and print result
solution = solve_with_bfs()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")