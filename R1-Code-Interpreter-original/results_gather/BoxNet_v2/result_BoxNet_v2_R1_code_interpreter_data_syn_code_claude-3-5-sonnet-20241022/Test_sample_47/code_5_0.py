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

def get_next_states(current_state):
    next_states = []
    occupied = set(current_state.values())
    
    # Try to move boxes that are not at their goals
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
            
        curr_pos = current_state[box]
        # Sort possible moves by distance to goal
        moves = sorted(adjacency[curr_pos],
                      key=lambda pos: manhattan_distance(pos, goal_state[box]))
        
        for next_pos in moves:
            if next_pos not in occupied:
                new_state = current_state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def solve():
    start_state = initial_state
    visited = {state_to_string(start_state)}
    
    # Priority queue entries: (heuristic, move_count, state_id, path)
    queue = [(0, 0, 0, [start_state])]
    state_counter = 1
    
    while queue:
        _, moves, _, path = heapq.heappop(queue)
        current_state = path[-1]
        
        # Check if we reached the goal
        if current_state == goal_state:
            return path
            
        # Generate and evaluate next states
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                
                # Calculate heuristic (sum of distances to goals)
                h = sum(manhattan_distance(next_state[box], goal_state[box])
                       for box in next_state)
                
                heapq.heappush(queue, (h + moves + 1, moves + 1, state_counter, new_path))
                state_counter += 1
    
    return None

# Solve and print result
solution = solve()
if solution:
    # Format the solution
    result = json.dumps(solution, separators=(',', ':'))
    print("<<<" + result + ">>>")
else:
    print("No solution found")