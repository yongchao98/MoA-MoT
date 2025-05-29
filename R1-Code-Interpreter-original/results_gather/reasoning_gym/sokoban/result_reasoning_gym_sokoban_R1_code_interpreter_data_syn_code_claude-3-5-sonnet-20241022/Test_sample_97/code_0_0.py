from collections import deque
import copy

def is_valid_pos(grid, row, col):
    return 0 <= row < len(grid) and 0 <= col < len(grid[0])

def get_next_states(grid, player_pos):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for direction, dy, dx in directions:
        new_y, new_x = player_pos[0] + dy, player_pos[1] + dx
        
        if not is_valid_pos(grid, new_y, new_x) or grid[new_y][new_x] == '+':
            continue
            
        new_grid = [list(row) for row in grid]
        new_player_pos = (new_y, new_x)
        
        # If moving to empty space or goal
        if grid[new_y][new_x] in '-X':
            if grid[player_pos[0]][player_pos[1]] == '*':
                new_grid[player_pos[0]][player_pos[1]] = 'X'
            else:
                new_grid[player_pos[0]][player_pos[1]] = '-'
            
            if grid[new_y][new_x] == 'X':
                new_grid[new_y][new_x] = '%'
            else:
                new_grid[new_y][new_x] = '*'
            next_states.append((direction, new_grid, new_player_pos))
            
        # If moving box
        elif grid[new_y][new_x] in '@$':
            push_y, push_x = new_y + dy, new_x + dx
            if (is_valid_pos(grid, push_y, push_x) and 
                grid[push_y][push_x] in '-X'):
                
                if grid[player_pos[0]][player_pos[1]] == '*':
                    new_grid[player_pos[0]][player_pos[1]] = 'X'
                else:
                    new_grid[player_pos[0]][player_pos[1]] = '-'
                
                if grid[new_y][new_x] == '$':
                    new_grid[new_y][new_x] = '%'
                else:
                    new_grid[new_y][new_x] = '*'
                
                if grid[push_y][push_x] == 'X':
                    new_grid[push_y][push_x] = '$'
                else:
                    new_grid[push_y][push_x] = '@'
                    
                next_states.append((direction, new_grid, new_player_pos))
    
    return next_states

def is_solved(grid):
    for row in grid:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def grid_to_string(grid):
    return ''.join(''.join(row) for row in grid)

def solve_sokoban(initial_grid):
    # Find player position
    player_pos = None
    for i in range(len(initial_grid)):
        for j in range(len(initial_grid[0])):
            if initial_grid[i][j] in '*%':
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([([], initial_grid, player_pos)])
    visited = {grid_to_string(initial_grid)}
    
    while queue:
        path, current_grid, current_pos = queue.popleft()
        
        if is_solved(current_grid):
            return ''.join(path)
        
        for direction, next_grid, next_pos in get_next_states(current_grid, current_pos):
            grid_str = grid_to_string(next_grid)
            if grid_str not in visited:
                visited.add(grid_str)
                queue.append((path + [direction], next_grid, next_pos))
    
    return None

# Parse input
puzzle = [
    "+++++++++++",
    "+---------+",
    "+-@XX--@--+",
    "+-$-@X---+",
    "+--X-----+",
    "+-$@+--+-+",
    "+---X+---+",
    "+$@@X--@-+",
    "+*----+X-+",
    "+++++++++++"
]

# Convert to list of lists for easier manipulation
grid = [list(row) for row in puzzle]
solution = solve_sokoban(grid)
print(solution)