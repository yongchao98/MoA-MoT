from collections import deque
import copy

def is_valid(grid, x, y):
    return 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != '+'

def get_next_states(grid, player_pos):
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for move, dx, dy in moves:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if is_valid(grid, new_x, new_y):
            new_grid = [row[:] for row in grid]
            if new_grid[new_x][new_y] in ['$', '@']:
                box_x, box_y = new_x + dx, new_y + dy
                if is_valid(grid, box_x, box_y) and grid[box_x][box_y] not in ['$', '@']:
                    new_grid[player_pos[0]][player_pos[1]] = '-'
                    new_grid[new_x][new_y] = '*'
                    new_grid[box_x][box_y] = '$'
                    next_states.append((move, new_grid, (new_x, new_y)))
            elif new_grid[new_x][new_y] == '-':
                new_grid[player_pos[0]][player_pos[1]] = '-'
                new_grid[new_x][new_y] = '*'
                next_states.append((move, new_grid, (new_x, new_y)))
    
    return next_states

def is_goal_state(grid):
    goals = [(1, 1), (4, 3)]
    boxes = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == '$':
                boxes.append((i, j))
    return all(box in boxes for box in goals)

def solve_sokoban(initial_grid, player_pos):
    queue = deque([([], initial_grid, player_pos)])
    visited = set()
    
    while queue:
        moves, grid, pos = queue.popleft()
        grid_tuple = tuple(map(tuple, grid))
        
        if (grid_tuple, pos) in visited:
            continue
            
        visited.add((grid_tuple, pos))
        
        if is_goal_state(grid):
            return ''.join(moves)
            
        for move, new_grid, new_pos in get_next_states(grid, pos):
            queue.append((moves + [move], new_grid, new_pos))
    
    return None

# Initial grid
grid = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', 'X', '-', '*', '-', '-', '-', '+'],
    ['+', '$', '@', '-', '-', '+', '$', '-', '+'],
    ['+', '-', '-', '-', '-', '+', '$', '-', '+'],
    ['+', '-', '-', '-', 'X', '@', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

player_pos = (1, 4)  # Starting position of the player
solution = solve_sokoban(grid, player_pos)
print(solution)