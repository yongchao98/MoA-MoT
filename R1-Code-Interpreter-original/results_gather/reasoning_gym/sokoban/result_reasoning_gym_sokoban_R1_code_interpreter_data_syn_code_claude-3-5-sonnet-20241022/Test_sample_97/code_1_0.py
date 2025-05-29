def get_player_pos(grid):
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in ['*', '%']:
                return (i, j)
    return None

def make_move(grid, player_pos, direction):
    y, x = player_pos
    dy, dx = direction
    new_y, new_x = y + dy, x + dx
    
    # Check if move is within bounds
    if not (0 <= new_y < len(grid) and 0 <= new_x < len(grid[0])):
        return None, None
    
    new_grid = [row[:] for row in grid]
    
    # Moving to wall
    if grid[new_y][new_x] == '+':
        return None, None
        
    # Moving to empty space or goal
    if grid[new_y][new_x] in '-X':
        new_grid[y][x] = 'X' if grid[y][x] == '*' else '-'
        new_grid[new_y][new_x] = '%' if grid[new_y][new_x] == 'X' else '*'
        return new_grid, (new_y, new_x)
        
    # Moving box
    if grid[new_y][new_x] in '@$':
        push_y, push_x = new_y + dy, new_x + dx
        if (0 <= push_y < len(grid) and 0 <= push_x < len(grid[0]) and 
            grid[push_y][push_x] in '-X'):
            new_grid[y][x] = 'X' if grid[y][x] == '*' else '-'
            new_grid[new_y][new_x] = '%' if grid[new_y][new_x] == '$' else '*'
            new_grid[push_y][push_x] = '$' if grid[push_y][push_x] == 'X' else '@'
            return new_grid, (new_y, new_x)
            
    return None, None

def is_solved(grid):
    return not any(cell in ['@', 'X'] for row in grid for cell in row)

def solve_sokoban(grid):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    visited = set()
    queue = [(grid, get_player_pos(grid), [])]
    max_moves = 30
    
    while queue:
        current_grid, pos, path = queue.pop(0)
        
        if len(path) > max_moves:
            continue
            
        if is_solved(current_grid):
            return ''.join(path)
            
        grid_str = str(current_grid)
        if grid_str in visited:
            continue
        visited.add(grid_str)
        
        for dir_char, dir_delta in directions.items():
            new_grid, new_pos = make_move(current_grid, pos, dir_delta)
            if new_grid is not None:
                queue.append((new_grid, new_pos, path + [dir_char]))
                
    return None

# Initialize puzzle
puzzle = [
    list("+++++++++++"),
    list("+---------+"),
    list("+-@XX--@--+"),
    list("+-$-@X---+"),
    list("+--X-----+"),
    list("+-$@+--+-+"),
    list("+---X+---+"),
    list("+$@@X--@-+"),
    list("+*----+X-+"),
    list("+++++++++++")
]

solution = solve_sokoban(puzzle)
print(solution)