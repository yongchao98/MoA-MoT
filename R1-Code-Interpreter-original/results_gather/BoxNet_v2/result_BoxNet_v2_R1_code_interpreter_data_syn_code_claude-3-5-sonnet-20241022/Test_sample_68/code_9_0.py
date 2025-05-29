import json
from collections import deque

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

def get_next_states(state):
    next_states = []
    # First try to move boxes that aren't at their goals
    boxes_to_move = []
    for box in state:
        if state[box] != goal_state[box]:
            boxes_to_move.append(box)
    
    # Sort boxes by distance to goal
    boxes_to_move.sort(key=lambda box: 
        abs(int(state[box].split(',')[0][1]) - int(goal_state[box].split(',')[0][1])) +
        abs(int(state[box].split(',')[1]) - int(goal_state[box].split(',')[1]))
    )
    
    # Generate moves for each box
    for box in boxes_to_move:
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
                
                # If this move reaches the goal for this box, prioritize it
                if next_pos == goal_state[box]:
                    return [new_state]
    
    return next_states

def solve():
    queue = deque([(initial_state, [initial_state])])
    visited = {str(sorted(initial_state.items()))}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
        
        for next_state in get_next_states(current_state):
            state_str = str(sorted(next_state.items()))
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
                
                if next_state == goal_state:
                    return new_path
    
    return None

# Get solution and format output
solution = solve()
if solution:
    # Format solution as list of sorted dictionaries
    formatted_solution = []
    for state in solution:
        formatted_state = dict(sorted(state.items()))
        formatted_solution.append(formatted_state)
    
    # Convert to string with specific formatting
    result = str(formatted_solution).replace("'", '"')
    # Print each state on a new line for better readability
    print(result)
else:
    print("No solution found")