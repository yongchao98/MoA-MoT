import heapq
import json
from copy import deepcopy
from dataclasses import dataclass
from typing import Dict, List

@dataclass
class SearchNode:
    priority: int
    cost: int
    path: List[Dict]
    state: Dict
    
    def __lt__(self, other):
        return self.priority < other.priority

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Simplified adjacency for the relevant paths
adjacency = {
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C4,7": ["C4,6", "C4,8", "C3,7"], "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"], "C4,1": ["C4,2", "C3,1"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    # Add intermediate positions
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"],
    "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"]
}

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def heuristic(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += manhattan_distance(state[box], goal_state[box])
    return total

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box in state:
        if state[box] == goal_state[box]:
            continue
        
        current_pos = state[box]
        if current_pos in adjacency:  # Check if position is in adjacency dict
            for next_pos in adjacency[current_pos]:
                if next_pos not in occupied:
                    new_state = deepcopy(state)
                    new_state[box] = next_pos
                    moves.append(new_state)
    
    return moves

def find_path(max_depth=50):
    start_state = initial_state
    visited = set()
    initial_node = SearchNode(
        priority=heuristic(start_state),
        cost=0,
        path=[start_state],
        state=start_state
    )
    pq = [initial_node]
    heapq.heapify(pq)
    
    while pq:
        current_node = heapq.heappop(pq)
        
        if len(current_node.path) > max_depth:
            continue
            
        state_str = state_to_string(current_node.state)
        
        if current_node.state == goal_state:
            return current_node.path
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state in get_valid_moves(current_node.state):
            if state_to_string(next_state) not in visited:
                new_path = current_node.path + [next_state]
                new_cost = current_node.cost + 1
                new_priority = new_cost + heuristic(next_state)
                new_node = SearchNode(
                    priority=new_priority,
                    cost=new_cost,
                    path=new_path,
                    state=next_state
                )
                heapq.heappush(pq, new_node)
    
    return None

# Find and print the solution
solution = find_path()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")