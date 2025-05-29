from collections import deque
import copy

def is_valid(grid, x, y):
    return 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != '+'

def get_next_states(grid, player_pos):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    next_states = []
    px, py = player_pos
    
    for dx, dy, move in directions:
        new_x, new_y = px + dx, py + dy
        
        if not is_valid(grid, new_x, new_y):
            continue
            
        new_grid = [list(row) for row in grid]
        
        # If moving to empty space or goal
        if grid[new_x][new_y] in ['-', 'X']:
            new_grid[px][py] = '-' if grid[px][py] == '*' else 'X'
            new_grid[new_x][new_y] = '*' if grid[new_x][new_y] == '-' else '%'
            next_states.append((new_grid, (new_x, new_y), move))
            
        # If moving box
        elif grid[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if is_valid(grid, box_x, box_y) and grid[box_x][box_y] in ['-', 'X']:
                new_grid[px][py] = '-' if grid[px][py] == '*' else 'X'
                new_grid[new_x][new_y] = '*' if grid[new_x][new_y] == '@' else '%'
                new_grid[box_x][box_y] = '@' if grid[box_x][box_y] == '-' else '$'
                next_states.append((new_grid, (new_x, new_y), move))
    
    return next_states

def is_solved(grid):
    for row in grid:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def grid_to_string(grid):
    return '\n'.join(''.join(row) for row in grid)

def solve_sokoban(initial_grid):
    # Find player position
    player_pos = None
    for i in range(len(initial_grid)):
        for j in range(len(initial_grid[0])):
            if initial_grid[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([(initial_grid, player_pos, "")])
    visited = {grid_to_string(initial_grid)}
    
    while queue:
        current_grid, current_pos, path = queue.popleft()
        
        if is_solved(current_grid):
            return path
            
        for next_grid, next_pos, move in get_next_states(current_grid, current_pos):
            grid_str = grid_to_string(next_grid)
            if grid_str not in visited:
                visited.add(grid_str)
                queue.append((next_grid, next_pos, path + move))
    
    return None

# Initial grid
initial_grid = [
    list("+++++++++++"),
    list("++X-@--++"),
    list("+++X-@$-X+"),
    list("++-@-$--+"),
    list("+-X++@+"),
    list("+-@---+"),
    list("+-----X+"),
    list("+-@*@X$+"),
    list("+++---X+"),
    list("+++++++++++")
]

solution = solve_sokoban(initial_grid)
print(solution)