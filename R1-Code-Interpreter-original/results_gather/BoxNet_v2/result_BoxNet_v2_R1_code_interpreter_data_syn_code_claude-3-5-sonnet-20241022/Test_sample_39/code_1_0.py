import copy
from collections import deque

initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C3,6"]
}

def find_path(start, goal, occupied_cells):
    if start == goal:
        return []
    
    visited = set()
    queue = deque([(start, [])])
    visited.add(start)
    
    while queue:
        current, path = queue.popleft()
        
        for next_cell in adjacency[current]:
            if next_cell not in visited and next_cell not in occupied_cells:
                if next_cell == goal:
                    return path + [next_cell]
                visited.add(next_cell)
                queue.append((next_cell, path + [next_cell]))
    return None

def is_goal_reached(current_state):
    return all(current_state[box] == goal_state[box] for box in current_state)

def generate_plan():
    states = [initial_state]
    current_state = copy.deepcopy(initial_state)
    max_steps = 100  # Prevent infinite loops
    step_count = 0
    
    while not is_goal_reached(current_state) and step_count < max_steps:
        moved = False
        for box in current_state:
            if current_state[box] == goal_state[box]:
                continue
                
            occupied_cells = {pos for b, pos in current_state.items() if b != box}
            path = find_path(current_state[box], goal_state[box], occupied_cells)
            
            if path:
                new_state = copy.deepcopy(current_state)
                new_state[box] = path[0]
                current_state = new_state
                states.append(new_state)
                moved = True
                break
        
        if not moved:
            return None
        step_count += 1
    
    return states if is_goal_reached(current_state) else None

# Generate and print the plan
plan = generate_plan()
if plan:
    print(plan)
else:
    print("No solution found")