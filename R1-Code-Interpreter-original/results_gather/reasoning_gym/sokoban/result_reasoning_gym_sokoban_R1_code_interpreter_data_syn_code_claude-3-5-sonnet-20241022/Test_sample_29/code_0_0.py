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

def solve_sokoban(grid):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    player = get_player_pos(grid)
    boxes = get_boxes(grid)
    goals = get_goals(grid)
    
    queue = deque([(player, boxes, "")])
    seen = set()
    
    while queue:
        pos, current_boxes, path = queue.popleft()
        
        if is_solved(current_boxes, goals):
            return path
            
        state = (pos, tuple(sorted(current_boxes)))
        if state in seen:
            continue
        seen.add(state)
        
        for move, (dr, dc) in directions.items():
            new_r, new_c = pos[0] + dr, pos[1] + dc
            
            if not is_valid(grid, new_r, new_c):
                continue
                
            if (new_r, new_c) in current_boxes:
                box_new_r, box_new_c = new_r + dr, new_c + dc
                if not is_valid(grid, box_new_r, box_new_c) or (box_new_r, box_new_c) in current_boxes:
                    continue
                new_boxes = set(current_boxes)
                new_boxes.remove((new_r, new_c))
                new_boxes.add((box_new_r, box_new_c))
                queue.append(((new_r, new_c), new_boxes, path + move))
            else:
                queue.append(((new_r, new_c), current_boxes, path + move))
    
    return None

# Input puzzle
grid = [
    list("+++++++++"),
    list("+-@-XX--+"),
    list("+-@X-@@-+"),
    list("+-@X---X+"),
    list("+-$-----+"),
    list("+*------+"),
    list("+++++++++")
]

solution = solve_sokoban(grid)
print(solution)