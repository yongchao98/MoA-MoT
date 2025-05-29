from collections import deque
import json

def find_shortest_path(adjacency, start, end):
    if start == end:
        return [start]
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        node = path[-1]
        
        for neighbor in adjacency[node]:
            if neighbor == end:
                return path + [neighbor]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(path + [neighbor])
    return None

def is_valid_move(state1, state2, adjacency):
    moved_boxes = []
    for box in state1:
        if state1[box] != state2[box]:
            if state2[box] not in adjacency[state1[box]]:
                return False
            moved_boxes.append((state1[box], state2[box]))
    
    # Check if more than one box moved
    if len(moved_boxes) > 1:
        return False
    
    # Check if any box moved to a cell occupied by another box
    occupied_cells = [state2[box] for box in state2]
    if len(occupied_cells) != len(set(occupied_cells)):
        return False
    
    return True

def solve_planning_problem():
    initial_state = {
        "box1": "C1,4", "box2": "C2,8",
        "box3": "C2,2", "box4": "C2,1"
    }
    goal_state = {
        "box1": "C2,8", "box2": "C2,2",
        "box3": "C1,2", "box4": "C2,5"
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
    
    solution = [initial_state]
    current_state = dict(initial_state)
    
    # First move box3 to C1,2
    path = find_shortest_path(adjacency, current_state["box3"], "C1,2")
    for next_pos in path[1:]:
        new_state = dict(current_state)
        new_state["box3"] = next_pos
        if is_valid_move(current_state, new_state, adjacency):
            solution.append(new_state)
            current_state = new_state
    
    # Move box4 to C2,5
    path = find_shortest_path(adjacency, current_state["box4"], "C2,5")
    for next_pos in path[1:]:
        new_state = dict(current_state)
        new_state["box4"] = next_pos
        if is_valid_move(current_state, new_state, adjacency):
            solution.append(new_state)
            current_state = new_state
    
    # Move box1 to C2,8
    path = find_shortest_path(adjacency, current_state["box1"], "C2,8")
    for next_pos in path[1:]:
        new_state = dict(current_state)
        new_state["box1"] = next_pos
        if is_valid_move(current_state, new_state, adjacency):
            solution.append(new_state)
            current_state = new_state
            
    # Finally move box2 to C2,2
    path = find_shortest_path(adjacency, current_state["box2"], "C2,2")
    for next_pos in path[1:]:
        new_state = dict(current_state)
        new_state["box2"] = next_pos
        if is_valid_move(current_state, new_state, adjacency):
            solution.append(new_state)
            current_state = new_state
    
    # Verify final state matches goal state
    if current_state != goal_state:
        return None
    
    print(json.dumps(solution))

solve_planning_problem()