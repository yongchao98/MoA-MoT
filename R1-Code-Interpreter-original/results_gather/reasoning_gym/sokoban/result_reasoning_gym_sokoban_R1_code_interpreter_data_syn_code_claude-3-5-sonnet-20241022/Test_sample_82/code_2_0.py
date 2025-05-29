from collections import deque

def is_valid(grid, row, col):
    if 0 <= row < len(grid) and 0 <= col < len(grid[0]):
        return grid[row][col] != '+'
    return False

def get_initial_state(grid):
    player = None
    boxes = set()
    goals = set()
    
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == '@':  # Player
                player = (i, j)
            elif grid[i][j] == '*':  # Player on empty space
                player = (i, j)
            elif grid[i][j] == '%':  # Player on goal
                player = (i, j)
                goals.add((i, j))
            elif grid[i][j] == '$':  # Box on goal
                boxes.add((i, j))
                goals.add((i, j))
            elif grid[i][j] == 'X':  # Goal
                goals.add((i, j))
            elif grid[i][j] == '@':  # Box
                boxes.add((i, j))
    
    return player, boxes, goals

def is_deadlock(grid, box_pos, goals):
    # Simple corner deadlock detection
    row, col = box_pos
    if box_pos not in goals:
        # Check if box is in corner
        if (not is_valid(grid, row-1, col) and not is_valid(grid, row, col-1)) or \
           (not is_valid(grid, row-1, col) and not is_valid(grid, row, col+1)) or \
           (not is_valid(grid, row+1, col) and not is_valid(grid, row, col-1)) or \
           (not is_valid(grid, row+1, col) and not is_valid(grid, row, col+1)):
            return True
    return False

def solve_sokoban(grid):
    player, boxes, goals = get_initial_state(grid)
    if not player or not boxes or not goals:
        return None
    
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    queue = deque([(player, boxes, "")])
    visited = set()
    
    while queue:
        curr_player, curr_boxes, path = queue.popleft()
        
        if curr_boxes == goals:
            return path
            
        state = (curr_player, frozenset(curr_boxes))
        if state in visited:
            continue
        visited.add(state)
        
        for move, dy, dx in directions:
            new_player_pos = (curr_player[0] + dy, curr_player[1] + dx)
            
            if not is_valid(grid, new_player_pos[0], new_player_pos[1]):
                continue
                
            new_boxes = set(curr_boxes)
            
            if new_player_pos in curr_boxes:  # If we're pushing a box
                new_box_pos = (new_player_pos[0] + dy, new_player_pos[1] + dx)
                
                if not is_valid(grid, new_box_pos[0], new_box_pos[1]) or \
                   new_box_pos in new_boxes or \
                   is_deadlock(grid, new_box_pos, goals):
                    continue
                    
                new_boxes.remove(new_player_pos)
                new_boxes.add(new_box_pos)
            
            queue.append((new_player_pos, new_boxes, path + move))
    
    return None

# Initialize the puzzle
puzzle = [
    list("+++++++" ),
    list("++$++++" ),
    list("+X+++$+" ),
    list("+-+++++" ),
    list("+---+$+" ),
    list("+--@+$+" ),
    list("+--X+++" ),
    list("+@--+++" ),
    list("+*@-X++" ),
    list("+++++++" )
]

solution = solve_sokoban(puzzle)
print(solution if solution else "No solution found")