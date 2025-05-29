from collections import deque
import heapq

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_goals_and_boxes(grid):
    goals = []
    boxes = []
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in ['X', '$', '%']:
                goals.append((i, j))
            if grid[i][j] in ['@', '$']:
                boxes.append((i, j))
    return goals, boxes

def heuristic(grid):
    goals, boxes = get_goals_and_boxes(grid)
    if len(boxes) != len(goals):
        return float('inf')
    
    total_dist = 0
    used_goals = set()
    
    for box in boxes:
        min_dist = float('inf')
        best_goal = None
        for goal in goals:
            if goal not in used_goals:
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    best_goal = goal
        if best_goal:
            used_goals.add(best_goal)
            total_dist += min_dist
    
    return total_dist

def get_next_states(grid, player_pos):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for direction, dy, dx in directions:
        new_y, new_x = player_pos[0] + dy, player_pos[1] + dx
        
        if not (0 <= new_y < len(grid) and 0 <= new_x < len(grid[0])) or grid[new_y][new_x] == '+':
            continue
            
        new_grid = [row[:] for row in grid]
        new_player_pos = (new_y, new_x)
        
        if grid[new_y][new_x] in '-X':
            new_grid[player_pos[0]][player_pos[1]] = 'X' if grid[player_pos[0]][player_pos[1]] == '*' else '-'
            new_grid[new_y][new_x] = '%' if grid[new_y][new_x] == 'X' else '*'
            next_states.append((direction, new_grid, new_player_pos))
            
        elif grid[new_y][new_x] in '@$':
            push_y, push_x = new_y + dy, new_x + dx
            if (0 <= push_y < len(grid) and 0 <= push_x < len(grid[0]) and 
                grid[push_y][push_x] in '-X'):
                
                new_grid[player_pos[0]][player_pos[1]] = 'X' if grid[player_pos[0]][player_pos[1]] == '*' else '-'
                new_grid[new_y][new_x] = '%' if grid[new_y][new_x] == '$' else '*'
                new_grid[push_y][push_x] = '$' if grid[push_y][push_x] == 'X' else '@'
                
                next_states.append((direction, new_grid, new_player_pos))
    
    return next_states

def is_solved(grid):
    return not any(cell in ['@', 'X'] for row in grid for cell in row)

def solve_sokoban(grid):
    player_pos = None
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    initial_h = heuristic(grid)
    pq = [(initial_h, 0, [], grid, player_pos)]
    visited = {str(grid)}
    max_steps = 50  # Limit search depth
    
    while pq:
        h, steps, path, current_grid, current_pos = heapq.heappop(pq)
        
        if is_solved(current_grid):
            return ''.join(path)
            
        if steps >= max_steps:
            continue
            
        for direction, next_grid, next_pos in get_next_states(current_grid, current_pos):
            grid_str = str(next_grid)
            if grid_str not in visited:
                visited.add(grid_str)
                next_h = heuristic(next_grid)
                heapq.heappush(pq, (next_h + steps + 1, steps + 1, path + [direction], next_grid, next_pos))
    
    return None

# Initialize puzzle
puzzle = [
    list("+++++++++++"),
    list("+---------+"),
    list("+-@XX--@--+"),
    list("+-$-@X---+"),
    list("+--X-----+"),
    list("+-$@+--+-+"),
    list("+---X+---+"),
    list("+$@@X--@-+"),
    list("+*----+X-+"),
    list("+++++++++++")
]

solution = solve_sokoban(puzzle)
print(solution)