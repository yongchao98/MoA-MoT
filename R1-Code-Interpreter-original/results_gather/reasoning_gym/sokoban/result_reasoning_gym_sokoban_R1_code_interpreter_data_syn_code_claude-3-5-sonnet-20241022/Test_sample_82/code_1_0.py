from collections import deque

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

def get_state_string(grid, player_pos, boxes):
    return str(player_pos) + str(sorted(list(boxes)))

def solve_sokoban(initial_grid):
    # Clean up the input grid
    grid = []
    for row in initial_grid:
        grid.append([c for c in row if c != ' '])
    
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    # Get initial state
    players = get_player_positions(grid)
    if not players:
        return None
    player = players[0]
    
    boxes = get_boxes(grid)
    goals = get_goals(grid)
    
    queue = deque([(player, boxes, "")])
    visited = set()
    
    while queue:
        current_player, current_boxes, path = queue.popleft()
        
        if is_solved(current_boxes, goals):
            return path
        
        state = get_state_string(grid, current_player, current_boxes)
        if state in visited:
            continue
        visited.add(state)
        
        for direction, (dy, dx) in directions.items():
            new_player = (current_player[0] + dy, current_player[1] + dx)
            
            if not is_valid(grid, new_player[0], new_player[1]):
                continue
                
            new_boxes = set(current_boxes)
            
            # If there's a box in the new position
            if new_player in new_boxes:
                box_new_pos = (new_player[0] + dy, new_player[1] + dx)
                
                if not is_valid(grid, box_new_pos[0], box_new_pos[1]) or box_new_pos in new_boxes:
                    continue
                    
                new_boxes.remove(new_player)
                new_boxes.add(box_new_pos)
            
            queue.append((new_player, new_boxes, path + direction))
    
    return None

# Initialize the puzzle
puzzle = [
    "+++++++" ,
    "++$++++" ,
    "+X+++$+" ,
    "+-+++++" ,
    "+---+$+" ,
    "+--@+$+" ,
    "+--X+++" ,
    "+@--+++" ,
    "+*@-X++" ,
    "+++++++"
]

solution = solve_sokoban(puzzle)
print(solution if solution else "No solution found")