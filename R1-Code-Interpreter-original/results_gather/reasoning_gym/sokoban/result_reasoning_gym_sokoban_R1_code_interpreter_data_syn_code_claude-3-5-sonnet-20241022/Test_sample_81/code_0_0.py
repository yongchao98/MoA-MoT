from collections import deque
import copy

def is_valid(grid, row, col):
    return 0 <= row < len(grid) and 0 <= col < len(grid[0]) and grid[row][col] != '+'

def get_next_pos(row, col, direction):
    if direction == 'U': return row-1, col
    if direction == 'D': return row+1, col
    if direction == 'L': return row, col-1
    if direction == 'R': return row, col+1

def convert_to_grid(puzzle):
    return [list(row.strip()) for row in puzzle.strip().split('\n')]

def find_player(grid):
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['*', '%']:
                return i, j
    return None

def find_boxes(grid):
    boxes = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['@', '$']:
                boxes.add((i, j))
    return boxes

def find_goals(grid):
    goals = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return goals

def solve_sokoban(puzzle):
    grid = convert_to_grid(puzzle)
    player = find_player(grid)
    boxes = find_boxes(grid)
    goals = find_goals(grid)
    
    queue = deque([(player, boxes, "")])
    visited = set()
    
    while queue:
        curr_player, curr_boxes, path = queue.popleft()
        
        if curr_boxes == goals:
            return path
            
        state = (curr_player, tuple(sorted(curr_boxes)))
        if state in visited:
            continue
        visited.add(state)
        
        for direction in ['U', 'D', 'L', 'R']:
            next_player = get_next_pos(*curr_player, direction)
            
            if not is_valid(grid, *next_player):
                continue
                
            if next_player in curr_boxes:
                next_box = get_next_pos(*next_player, direction)
                if not is_valid(grid, *next_box) or next_box in curr_boxes:
                    continue
                    
                new_boxes = set(curr_boxes)
                new_boxes.remove(next_player)
                new_boxes.add(next_box)
                queue.append((next_player, new_boxes, path + direction))
            else:
                queue.append((next_player, curr_boxes, path + direction))
    
    return None

puzzle = """+ + + + + +
+ X - - - +
+ - @ - $ +
+ - $ + + +
+ - @ X - +
+ - * $ - +
+ - - - - +
+ + + + + +"""

solution = solve_sokoban(puzzle)
print(solution)