import json
from collections import deque
import heapq

def heuristic(cell1, cell2):
    x1, y1 = map(int, cell1.replace('C', '').split(','))
    x2, y2 = map(int, cell2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path(came_from, current):
    path = [current]
    while current in came_from:
        current = came_from[current]
        path.append(current)
    return path[::-1]

def a_star(start, goal, adjacency, obstacles):
    frontier = [(0, start)]
    came_from = {}
    cost_so_far = {start: 0}
    
    while frontier:
        _, current = heapq.heappop(frontier)
        
        if current == goal:
            return get_path(came_from, current)
            
        for next_cell in adjacency[current]:
            if next_cell in obstacles:
                continue
                
            new_cost = cost_so_far[current] + 1
            if next_cell not in cost_so_far or new_cost < cost_so_far[next_cell]:
                cost_so_far[next_cell] = new_cost
                priority = new_cost + heuristic(next_cell, goal)
                heapq.heappush(frontier, (priority, next_cell))
                came_from[next_cell] = current
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    unvisited_goals = set(goals)
    path = [start]
    
    while unvisited_goals:
        min_dist = float('inf')
        next_goal = None
        best_path = None
        
        for goal in unvisited_goals:
            temp_path = a_star(current, goal, adjacency, obstacles)
            if temp_path and len(temp_path) < min_dist:
                min_dist = len(temp_path)
                next_goal = goal
                best_path = temp_path
        
        if not next_goal:
            return None
            
        path.extend(best_path[1:])
        current = next_goal
        unvisited_goals.remove(next_goal)
    
    return path

# Problem setup
obstacles = ['C5,1', 'C1,1', 'C3,1', 'C3,3', 'C2,2', 'C4,5', 'C5,4', 'C2,4']
goals = ['C3,5', 'C2,5', 'C5,3', 'C1,4', 'C1,2', 'C4,2']
start = 'C4,4'

# Complete adjacency dictionary
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

# Find solution
solution = find_complete_path(start, goals, adjacency, obstacles)
print(json.dumps(solution))