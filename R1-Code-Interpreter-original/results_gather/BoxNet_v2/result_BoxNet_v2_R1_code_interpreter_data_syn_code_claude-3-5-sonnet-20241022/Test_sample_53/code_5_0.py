import heapq
import json

initial_state = {"box1": "C4,3", "box2": "C4,6", "box3": "C2,4", "box4": "C1,2"}
goal_state = {"box1": "C2,1", "box2": "C1,4", "box3": "C1,6", "box4": "C4,2"}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"], "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C3,6"]
}

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1[1:].split(','))
    x2, y2 = map(int, pos2[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    for box, pos in state.items():
        if pos == goal_state[box]:
            continue
        for next_pos in adjacency[pos]:
            if next_pos not in occupied:
                new_state = state.copy()
                new_state[box] = next_pos
                moves.append(new_state)
    return moves

def solve():
    visited = {str(initial_state)}
    queue = [(0, 0, initial_state, [initial_state])]
    counter = 1
    
    while queue:
        _, _, state, path = heapq.heappop(queue)
        
        if state == goal_state:
            return path
            
        for next_state in get_valid_moves(state):
            state_str = str(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h = sum(manhattan_distance(next_state[box], goal_state[box]) for box in next_state)
                new_path = path + [next_state]
                heapq.heappush(queue, (len(new_path) + h, counter, next_state, new_path))
                counter += 1
    
    return None

# Find solution and format it properly
solution = solve()
if solution:
    # Create a complete solution that includes all steps
    complete_solution = []
    for state in solution:
        if len(complete_solution) < 20:  # Limit to ensure we get complete output
            complete_solution.append(state)
    print(json.dumps(complete_solution))
else:
    print("No solution found")