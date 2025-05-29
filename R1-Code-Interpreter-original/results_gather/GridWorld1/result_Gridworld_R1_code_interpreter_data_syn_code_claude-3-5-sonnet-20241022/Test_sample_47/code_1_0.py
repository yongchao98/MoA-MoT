from collections import deque
import json

def find_path(start, goals, obstacles, adjacency):
    def bfs(start, target, visited_set=None):
        if visited_set is None:
            visited_set = set()
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
                
            for next_pos in adjacency.get(current, []):
                if next_pos not in visited and next_pos not in obstacles:
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    def find_best_sequence(current, remaining_goals, path_so_far):
        if not remaining_goals:
            return path_so_far
            
        best_path = None
        for goal in remaining_goals:
            # Find path to this goal
            segment = bfs(current, goal)
            if segment:
                # Try this goal
                new_path = path_so_far + segment[1:]  # Exclude start position to avoid duplicates
                result = find_best_sequence(
                    goal,
                    [g for g in remaining_goals if g != goal],
                    new_path
                )
                if result and (best_path is None or len(result) < len(best_path)):
                    best_path = result
                    
        return best_path

    # Start the search
    initial_path = [start]
    final_path = find_best_sequence(start, goals, initial_path)
    
    return final_path

# Problem setup
start = "C1,5"
goals = ['C2,2', 'C6,1', 'C6,3', 'C5,3', 'C5,4', 'C3,5', 'C2,4']
obstacles = ['C1,6', 'C2,3', 'C6,5', 'C5,1', 'C5,6', 'C2,6', 'C6,4', 'C1,2', 'C5,2', 'C3,4']

# Define adjacency dictionary
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

# Find the path
path = find_path(start, goals, set(obstacles), adjacency)

# Format and print the result
if path:
    result = json.dumps(path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")