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
    
    queue = deque([(start, [])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in occupied_cells:
                if next_pos == goal:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def generate_plan():
    states = [initial_state]
    current_state = copy.deepcopy(initial_state)
    boxes_to_move = set(initial_state.keys())
    
    while boxes_to_move:
        moved = False
        for box in list(boxes_to_move):
            if current_state[box] == goal_state[box]:
                boxes_to_move.remove(box)
                continue
            
            # Get occupied cells excluding current box
            occupied_cells = {pos for b, pos in current_state.items() if b != box}
            
            # Try to find a path to the goal
            path = find_path(current_state[box], goal_state[box], occupied_cells)
            
            if path:
                # Move box one step along the path
                new_state = copy.deepcopy(current_state)
                new_state[box] = path[0]
                current_state = new_state
                states.append(new_state)
                moved = True
                break
        
        if not moved and boxes_to_move:
            # If we can't move any box directly towards its goal,
            # try moving boxes to make way
            for box in list(boxes_to_move):
                occupied_cells = {pos for b, pos in current_state.items() if b != box}
                for adj_pos in adjacency[current_state[box]]:
                    if adj_pos not in occupied_cells:
                        new_state = copy.deepcopy(current_state)
                        new_state[box] = adj_pos
                        current_state = new_state
                        states.append(new_state)
                        moved = True
                        break
                if moved:
                    break
            
            if not moved:
                return None
        
        if len(states) > 100:  # Prevent infinite loops
            return None
    
    return states

# Generate and print the plan
plan = generate_plan()
if plan:
    print(plan)
else:
    print("No solution found")