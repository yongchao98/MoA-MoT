import json
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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def move_box_to_goal(states, box, path, occupied):
    current = states[-1][box]
    target = goal_state[box]
    
    if current == target:
        return True
        
    if not path:
        # Find path to goal using simple BFS
        visited = {current}
        queue = [(current, [current])]
        while queue:
            pos, p = queue.pop(0)
            for next_pos in adjacency[pos]:
                if next_pos not in visited and next_pos not in occupied:
                    if next_pos == target:
                        path.extend(p + [next_pos])
                        break
                    visited.add(next_pos)
                    queue.append((next_pos, p + [next_pos]))
            if path:
                break
                
    if path:
        next_pos = path.pop(0)
        new_state = deepcopy(states[-1])
        new_state[box] = next_pos
        states.append(new_state)
        return True
        
    return False

def solve():
    states = [initial_state]
    paths = {"box1": [], "box2": [], "box3": [], "box4": []}
    
    while states[-1] != goal_state:
        moved = False
        current = states[-1]
        occupied = set(current.values())
        
        # Try to move each box that's not at its goal
        for box in ["box1", "box2", "box3", "box4"]:
            if current[box] != goal_state[box]:
                occupied.remove(current[box])
                if move_box_to_goal(states, box, paths[box], occupied):
                    moved = True
                    break
                occupied.add(current[box])
        
        if not moved:
            # If we can't move any box, we're stuck
            return None
    
    return states

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))