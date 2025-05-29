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

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) * 2
              for box in state if state[box] != goal_state[box])

def get_next_states(current_state):
    next_states = []
    occupied = set(current_state.values())
    
    # Prioritize boxes that are closer to their goals
    boxes = sorted(current_state.keys(),
                  key=lambda b: manhattan_distance(current_state[b], goal_state[b]))
    
    for box in boxes:
        if current_state[box] == goal_state[box]:
            continue
        
        curr_pos = current_state[box]
        # Prioritize moves towards the goal
        possible_moves = sorted(adjacency[curr_pos],
                              key=lambda pos: manhattan_distance(pos, goal_state[box]))
        
        for next_pos in possible_moves:
            if next_pos not in occupied:
                new_state = current_state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
                
                # If this move gets a box to its goal, prioritize it
                if next_pos == goal_state[box]:
                    return [new_state]
    
    return next_states

def solve():
    start = initial_state
    visited = {json.dumps(start, sort_keys=True)}
    queue = [(heuristic(start), 0, 0, [start])]
    
    while queue:
        _, cost, _, path = heapq.heappop(queue)
        current = path[-1]
        
        if current == goal_state:
            return path
        
        if cost > 30:  # Limit path length to avoid excessive solutions
            continue
            
        for next_state in get_next_states(current):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                new_cost = len(new_path) - 1
                priority = new_cost + heuristic(next_state)
                heapq.heappush(queue, (priority, new_cost, len(visited), new_path))
    
    return None

solution = solve()
if solution:
    # Format output with readable indentation
    result = json.dumps(solution)
    print("<<<" + result + ">>>")
else:
    print("No solution found")