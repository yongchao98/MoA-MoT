from collections import deque, defaultdict
import json
import copy

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
    
    # Check if more than one box moved from or to the same cell
    cells_from = [x[0] for x in moved_boxes]
    cells_to = [x[1] for x in moved_boxes]
    if len(cells_from) != len(set(cells_from)) or len(cells_to) != len(set(cells_to)):
        return False
    
    return True

def solve_planning_problem():
    # Initial and goal states
    initial_state = {
        "box1": "C1,4", "box2": "C2,8",
        "box3": "C2,2", "box4": "C2,1"
    }
    goal_state = {
        "box1": "C2,8", "box2": "C2,2",
        "box3": "C1,2", "box4": "C2,5"
    }
    
    # Adjacency information
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
    
    # Find paths for each box
    solution = [initial_state]
    current_state = dict(initial_state)
    
    # Order boxes by distance to goal
    boxes_order = []
    for box in initial_state:
        if initial_state[box] != goal_state[box]:
            path = find_shortest_path(adjacency, initial_state[box], goal_state[box])
            boxes_order.append((box, len(path) if path else float('inf')))
    
    boxes_order.sort(key=lambda x: x[1])
    
    # Move boxes one by one
    for box, _ in boxes_order:
        current_pos = current_state[box]
        target_pos = goal_state[box]
        
        if current_pos == target_pos:
            continue
            
        path = find_shortest_path(adjacency, current_pos, target_pos)
        
        for next_pos in path[1:]:
            new_state = dict(current_state)
            new_state[box] = next_pos
            
            if is_valid_move(current_state, new_state, adjacency):
                solution.append(new_state)
                current_state = new_state
    
    print(json.dumps(solution))

solve_planning_problem()