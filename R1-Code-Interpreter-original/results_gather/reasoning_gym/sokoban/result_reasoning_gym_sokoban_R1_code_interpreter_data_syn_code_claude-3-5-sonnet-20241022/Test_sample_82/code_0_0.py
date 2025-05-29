from collections import deque
import copy

def is_valid(grid, row, col):
    return 0 <= row < len(grid) and 0 <= col < len(grid[0]) and grid[row][col] != '+'

def get_player_positions(grid):
    players = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['*', '@', '%']:
                players.append((i, j))
    return players

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

def get_state_string(grid, player_pos):
    state = ""
    for row in grid:
        state += "".join(row)
    state += f",{player_pos}"
    return state

def solve_sokoban(initial_grid):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    grid = [list(row) for row in initial_grid]
    
    players = get_player_positions(grid)
    if not players:
        return None
    player = players[0]
    
    boxes = get_boxes(grid)
    goals = get_goals(grid)
    
    queue = deque([(grid, player, "", set())])
    visited = set()
    
    while queue:
        current_grid, current_player, path, current_boxes = queue.popleft()
        
        if is_solved(current_boxes, goals):
            return path
        
        state = get_state_string(current_grid, current_player)
        if state in visited:
            continue
        visited.add(state)
        
        for direction, (dy, dx) in directions.items():
            new_player = (current_player[0] + dy, current_player[1] + dx)
            
            if not is_valid(current_grid, new_player[0], new_player[1]):
                continue
                
            new_grid = [row[:] for row in current_grid]
            new_boxes = set(current_boxes)
            
            # If there's a box in the new position
            if (new_player[0], new_player[1]) in new_boxes:
                box_new_pos = (new_player[0] + dy, new_player[1] + dx)
                
                if not is_valid(current_grid, box_new_pos[0], box_new_pos[1]) or box_new_pos in new_boxes:
                    continue
                    
                new_boxes.remove(new_player)
                new_boxes.add(box_new_pos)
            
            queue.append((new_grid, new_player, path + direction, new_boxes))
    
    return None

# Initialize the puzzle
puzzle = [
    "+++++++" ,
    "++$++++" ,
    "+X+++$+" ,
    "+-++++" ,
    "+---+$+" ,
    "+--@+$+" ,
    "+--X+++" ,
    "+@--+++" ,
    "+*@-X++" ,
    "+++++++"
]

solution = solve_sokoban(puzzle)
print(solution if solution else "No solution found")