import json
from heapq import heappush, heappop

initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

def get_adjacent(pos):
    r, c = int(pos[1]), int(pos[3])
    adj = []
    if r > 1: adj.append(f"C{r-1},{c}")
    if r < 4: adj.append(f"C{r+1},{c}")
    if c > 1: adj.append(f"C{r},{c-1}")
    if c < 8: adj.append(f"C{r},{c+1}")
    return adj

def manhattan(pos1, pos2):
    r1, c1 = int(pos1[1]), int(pos1[3])
    r2, c2 = int(pos2[1]), int(pos2[3])
    return abs(r1 - r2) + abs(c1 - c2)

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def solve():
    visited = set()
    queue = [(0, 0, initial_state, [initial_state])]  # (priority, step, state, path)
    visited.add(state_to_string(initial_state))
    
    while queue:
        _, step, current, path = heappop(queue)
        
        if current == goal_state:
            return path
            
        if step > 20:  # Limit path length
            continue
            
        # Try moving each box
        for box in current:
            if current[box] == goal_state[box]:
                continue
                
            pos = current[box]
            for new_pos in get_adjacent(pos):
                if new_pos not in current.values():
                    next_state = current.copy()
                    next_state[box] = new_pos
                    state_str = state_to_string(next_state)
                    
                    if state_str not in visited:
                        visited.add(state_str)
                        h_score = sum(manhattan(next_state[b], goal_state[b]) for b in next_state)
                        priority = step + 1 + h_score
                        heappush(queue, (priority, step + 1, next_state, path + [next_state]))
    
    return None

# Find solution
solution = solve()
if solution:
    # Format solution
    result = []
    for state in solution:
        result.append(dict(sorted(state.items())))
    print(json.dumps(result))
else:
    print("No solution found")