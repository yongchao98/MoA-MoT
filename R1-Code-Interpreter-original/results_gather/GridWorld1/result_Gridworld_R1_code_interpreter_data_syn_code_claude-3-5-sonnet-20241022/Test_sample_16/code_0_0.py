import json
from typing import List, Set, Dict

def find_path(
    start: str,
    goals: List[str],
    obstacles: List[str],
    adjacency: Dict[str, List[str]]
) -> List[str]:
    def is_valid_move(pos: str) -> bool:
        return pos not in obstacles

    def dfs(
        current: str,
        unvisited_goals: Set[str],
        path: List[str],
        visited: Set[str]
    ) -> List[str]:
        # If we've visited all goals, we're done
        if not unvisited_goals:
            return path
        
        # Try each possible next position
        for next_pos in adjacency[current]:
            if not is_valid_move(next_pos):
                continue
                
            # Skip if we've been here before, unless it's a goal we still need
            if next_pos in visited and next_pos not in unvisited_goals:
                continue
                
            new_path = path + [next_pos]
            new_visited = visited | {next_pos}
            new_unvisited = unvisited_goals - {next_pos}
            
            result = dfs(next_pos, new_unvisited, new_path, new_visited)
            if result:
                return result
        
        return None

    # Initialize the search
    initial_path = [start]
    initial_visited = {start}
    goals_set = set(goals)
    initial_unvisited = goals_set - {start}
    
    # Find and return the path
    result = dfs(start, initial_unvisited, initial_path, initial_visited)
    return result

# Problem inputs
start_pos = "C3,4"
goals = ['C3,1', 'C5,2', 'C4,4', 'C1,2', 'C5,1', 'C4,3']
obstacles = ['C2,5', 'C3,5', 'C3,2', 'C2,2', 'C4,2', 'C2,4', 'C5,5', 'C4,5']

# Load adjacency data
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

# Find the path
path = find_path(start_pos, goals, obstacles, adjacency)

# Output the result
print(json.dumps(path))