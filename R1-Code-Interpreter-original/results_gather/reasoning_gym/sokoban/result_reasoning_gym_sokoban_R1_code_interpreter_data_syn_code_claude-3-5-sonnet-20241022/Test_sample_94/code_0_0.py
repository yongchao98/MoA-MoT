from collections import deque

def is_valid(grid, row, col):
    return 0 <= row < len(grid) and 0 <= col < len(grid[0]) and grid[row][col] != '+'

def get_next_states(state, grid):
    player_pos, boxes = state
    moves = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    next_states = []
    
    for dr, dc, move in moves:
        new_r, new_c = player_pos[0] + dr, player_pos[1] + dc
        
        if not is_valid(grid, new_r, new_c):
            continue
            
        new_boxes = boxes.copy()
        if (new_r, new_c) in boxes:
            push_r, push_c = new_r + dr, new_c + dc
            if not is_valid(grid, push_r, push_c) or (push_r, push_c) in boxes:
                continue
            new_boxes.remove((new_r, new_c))
            new_boxes.add((push_r, push_c))
        
        next_states.append(((new_r, new_c), frozenset(new_boxes), move))
    
    return next_states

def solve_sokoban(grid):
    rows, cols = len(grid), len(grid[0])
    goals = set()
    boxes = set()
    player_pos = None
    
    # Parse initial state
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ['X', '*', '$']:
                goals.add((i, j))
            if grid[i][j] in ['@', '*']:
                player_pos = (i, j)
            if grid[i][j] in ['@', '$']:
                boxes.add((i, j))
    
    # BFS
    start_state = (player_pos, frozenset(boxes))
    queue = deque([(start_state, "")])
    visited = {start_state}
    
    while queue:
        state, path = queue.popleft()
        
        # Check if all boxes are on goals
        if all(box in goals for box in state[1]):
            return path
        
        for next_state in get_next_states(state, grid):
            new_pos, new_boxes, move = next_state
            new_state = (new_pos, new_boxes)
            
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move))
    
    return None

# Initialize the grid
grid = [
    list("+++++++++++"),
    list("++-----+-+"),
    list("++X-@-@-+"),
    list("++X$@---+"),
    list("++$*$-X-+"),
    list("+++++++++++")
]

solution = solve_sokoban(grid)
print(solution)