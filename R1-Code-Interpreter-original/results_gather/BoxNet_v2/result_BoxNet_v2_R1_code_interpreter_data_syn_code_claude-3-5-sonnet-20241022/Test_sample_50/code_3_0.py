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

def get_next_states(state):
    states = []
    for box in state:
        if state[box] == goal_state[box]:
            continue
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in state.values():
                new_state = copy.deepcopy(state)
                new_state[box] = next_pos
                states.append(new_state)
    return states

def solve():
    start = initial_state
    queue = [(start, [start])]
    visited = {json.dumps(start, sort_keys=True)}
    
    while queue:
        current_state, path = queue.pop(0)
        
        # Check if we've reached the goal
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
            
        # Try moving one box at a time
        for next_state in get_next_states(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
                
                # Early stopping if we find a state where more boxes are in their goals
                current_goals = sum(1 for box in current_state if current_state[box] == goal_state[box])
                next_goals = sum(1 for box in next_state if next_state[box] == goal_state[box])
                if next_goals > current_goals:
                    queue.insert(0, (next_state, new_path))
    
    return None

# Solve and print result
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")