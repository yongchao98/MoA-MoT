import json
from collections import deque
from copy import deepcopy

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
    occupied = set(state.values())
    
    for box in state:
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                states.append(new_state)
    return states

def bfs():
    start = initial_state
    queue = deque([(start, [start])])
    visited = {json.dumps(start, sort_keys=True)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we've reached the goal
        if all(current_state[box] == goal_state[box] for box in goal_state):
            return path
            
        # Get all possible next states
        for next_state in get_next_states(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                
                # Calculate how many boxes are in their goal positions
                goal_boxes = sum(1 for box in next_state if next_state[box] == goal_state[box])
                current_goal_boxes = sum(1 for box in current_state if current_state[box] == goal_state[box])
                
                # If this move gets more boxes to their goals, prioritize it
                if goal_boxes > current_goal_boxes:
                    queue.appendleft((next_state, path + [next_state]))
                else:
                    queue.append((next_state, path + [next_state]))
                
        # Limit search space
        if len(visited) > 5000:
            return None

def solve():
    # Try to solve with initial approach
    solution = bfs()
    
    if solution:
        # Remove consecutive duplicates
        simplified = [solution[0]]
        for state in solution[1:]:
            if state != simplified[-1]:
                simplified.append(state)
        return simplified
    
    return None

# Solve and print result
solution = solve()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")