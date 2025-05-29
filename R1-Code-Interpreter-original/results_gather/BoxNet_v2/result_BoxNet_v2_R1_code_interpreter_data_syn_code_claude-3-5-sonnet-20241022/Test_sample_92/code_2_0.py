from queue import PriorityQueue
import json
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

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def heuristic(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += manhattan_distance(state[box], goal_state[box])
    return total

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    # Sort boxes to prioritize those furthest from their goals
    boxes = sorted(current_state.keys(), 
                  key=lambda b: -manhattan_distance(current_state[b], goal_state[b]))
    
    for box in boxes:
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        # Sort adjacent positions by distance to goal
        adj_positions = sorted(adjacency[current_pos],
                             key=lambda p: manhattan_distance(p, goal_state[box]))
        
        for next_pos in adj_positions:
            if next_pos not in occupied_cells:
                new_state = copy.deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    start_state = initial_state
    frontier = PriorityQueue()
    counter = 0
    frontier.put((0, counter, [start_state]))
    visited = {state_to_string(start_state)}
    
    while not frontier.empty():
        _, _, path = frontier.get()
        current_state = path[-1]
        
        if current_state == goal_state:
            return path
            
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                priority = len(new_path) + heuristic(next_state)
                counter += 1
                frontier.put((priority, counter, new_path))
    
    return None

solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")