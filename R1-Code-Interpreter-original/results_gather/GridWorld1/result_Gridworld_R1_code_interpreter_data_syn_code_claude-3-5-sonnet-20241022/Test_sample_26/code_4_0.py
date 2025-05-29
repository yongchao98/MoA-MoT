from collections import deque
import json

def find_shortest_path(adjacency, start, end, obstacles):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def get_distance(path):
    return len(path) if path else float('inf')

def find_complete_path(adjacency, start, goals, obstacles):
    current_pos = start
    remaining_goals = set(goals)
    final_path = [start]
    visited_goals = set()
    
    while remaining_goals:
        # Find all possible next goals and their paths
        possible_moves = []
        for goal in remaining_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path:
                # Calculate additional goals that would be visited along this path
                bonus_goals = len([p for p in path if p in remaining_goals])
                score = len(path) - bonus_goals * 0.5  # Favor paths that visit multiple goals
                possible_moves.append((score, path, goal))
        
        if not possible_moves:
            return None
        
        # Choose the most efficient next move
        possible_moves.sort(key=lambda x: x[0])  # Sort by score
        _, next_path, next_goal = possible_moves[0]
        
        # Add path to final path (excluding start position)
        final_path.extend(next_path[1:])
        
        # Update current position and remove all goals that were visited
        current_pos = next_goal
        for pos in next_path:
            if pos in remaining_goals:
                remaining_goals.remove(pos)
                visited_goals.add(pos)
    
    # Optimize the final path by removing any unnecessary moves
    optimized_path = [final_path[0]]
    for i in range(1, len(final_path)):
        current = final_path[i]
        if (current in goals or 
            i == len(final_path) - 1 or 
            current != optimized_path[-1]):
            optimized_path.append(current)
    
    return optimized_path

# Problem setup
start = "C2,6"
goals = ['C6,3', 'C4,2', 'C2,5', 'C4,3', 'C4,5', 'C1,1', 'C3,3']
obstacles = ['C5,6', 'C4,6', 'C1,4', 'C1,3', 'C6,2', 'C3,4', 'C6,1', 'C3,2', 'C4,1', 'C6,6']

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
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
    "C4,6": ["C4,5", "C3,6", "C5,6"],
    "C5,1": ["C5,2", "C4,1", "C6,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
    "C5,6": ["C5,5", "C4,6", "C6,6"],
    "C6,1": ["C6,2", "C5,1"],
    "C6,2": ["C6,1", "C6,3", "C5,2"],
    "C6,3": ["C6,2", "C6,4", "C5,3"],
    "C6,4": ["C6,3", "C6,5", "C5,4"],
    "C6,5": ["C6,4", "C6,6", "C5,5"],
    "C6,6": ["C6,5", "C5,6"]
}

# Find the solution path
solution_path = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if solution_path:
    print(f"<<<{json.dumps(solution_path)}>>>")
else:
    print("No valid path found")