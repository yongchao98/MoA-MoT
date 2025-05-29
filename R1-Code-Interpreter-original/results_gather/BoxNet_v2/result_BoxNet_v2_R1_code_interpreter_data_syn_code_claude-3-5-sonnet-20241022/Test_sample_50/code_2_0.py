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

def solve_sequential():
    current_state = copy.deepcopy(initial_state)
    solution = [copy.deepcopy(current_state)]
    
    # Move box1 to its goal
    path = move_box_to_goal(current_state, "box1")
    solution.extend(path[1:])
    current_state = path[-1]
    
    # Move box2 to its goal
    path = move_box_to_goal(current_state, "box2")
    solution.extend(path[1:])
    current_state = path[-1]
    
    # Move box3 to its goal
    path = move_box_to_goal(current_state, "box3")
    solution.extend(path[1:])
    current_state = path[-1]
    
    # Move box4 to its goal
    path = move_box_to_goal(current_state, "box4")
    solution.extend(path[1:])
    
    return solution

def move_box_to_goal(state, box_name):
    if state[box_name] == goal_state[box_name]:
        return [state]
        
    visited = {json.dumps(state, sort_keys=True)}
    queue = [(state, [state])]
    
    while queue:
        current_state, path = queue.pop(0)
        
        if current_state[box_name] == goal_state[box_name]:
            return path
            
        current_pos = current_state[box_name]
        for next_pos in adjacency[current_pos]:
            if next_pos not in current_state.values():
                next_state = copy.deepcopy(current_state)
                next_state[box_name] = next_pos
                state_str = json.dumps(next_state, sort_keys=True)
                
                if state_str not in visited:
                    visited.add(state_str)
                    queue.append((next_state, path + [next_state]))
    
    return None

# Solve and print result
solution = solve_sequential()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")