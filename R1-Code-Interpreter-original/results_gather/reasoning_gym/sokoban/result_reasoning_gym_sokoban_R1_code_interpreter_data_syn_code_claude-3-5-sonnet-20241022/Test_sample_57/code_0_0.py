from collections import deque
import copy

def is_valid(grid, r, c):
    return 0 <= r < len(grid) and 0 <= c < len(grid[0]) and grid[r][c] != '+'

def get_player_pos(grid):
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['*', '%']:
                return (i, j)
    return None

def get_boxes(grid):
    boxes = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def get_goals(grid):
    goals = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return goals

def is_solved(boxes, goals):
    return boxes == goals

def get_state(player, boxes):
    return (player, frozenset(boxes))

def solve_sokoban(grid):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    player = get_player_pos(grid)
    boxes = get_boxes(grid)
    goals = get_goals(grid)
    
    queue = deque([(player, boxes, "")])
    visited = {get_state(player, boxes)}
    
    while queue:
        player, boxes, path = queue.popleft()
        
        if is_solved(boxes, goals):
            return path
            
        for move, (dr, dc) in directions.items():
            new_r, new_c = player[0] + dr, player[1] + dc
            
            if not is_valid(grid, new_r, new_c):
                continue
                
            if (new_r, new_c) in boxes:
                box_r, box_c = new_r + dr, new_c + dc
                if not is_valid(grid, box_r, box_c) or (box_r, box_c) in boxes:
                    continue
                    
                new_boxes = set(boxes)
                new_boxes.remove((new_r, new_c))
                new_boxes.add((box_r, box_c))
                
                new_state = get_state((new_r, new_c), new_boxes)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(((new_r, new_c), new_boxes, path + move))
            else:
                new_state = get_state((new_r, new_c), boxes)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(((new_r, new_c), boxes, path + move))
    
    return None

# Initialize the grid
grid = [
    list("++++++"),
    list("+---*+"),
    list("+-X@-+"),
    list("+--@X+"),
    list("+--@-+"),
    list("+-@X++"),
    list("+$X+++"),
    list("++++++"),
    list("+$++++"),
    list("++++++")
]

solution = solve_sokoban(grid)
print(solution)