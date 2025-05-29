from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(grid, players):
    rows, cols = len(grid), len(grid[0])
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for player_idx, (px, py) in enumerate(players):
        for move, dx, dy in moves:
            new_x, new_y = px + dx, py + dy
            
            if not is_valid(new_x, new_y, rows, cols) or grid[new_x][new_y] == '+':
                continue
                
            new_grid = [row[:] for row in grid]
            new_players = players[:]
            
            # If moving to empty space or goal
            if grid[new_x][new_y] in ['-', 'X']:
                new_grid[px][py] = '-' if grid[px][py] not in ['X', '$'] else 'X'
                new_grid[new_x][new_y] = '*' if grid[new_x][new_y] == '-' else '%'
                new_players[player_idx] = (new_x, new_y)
                next_states.append((move, new_grid, new_players))
                
            # If moving a box
            elif grid[new_x][new_y] in ['@', '$']:
                box_x, box_y = new_x + dx, new_y + dy
                if (is_valid(box_x, box_y, rows, cols) and 
                    grid[box_x][box_y] in ['-', 'X']):
                    new_grid[px][py] = '-' if grid[px][py] not in ['X', '$'] else 'X'
                    new_grid[new_x][new_y] = '*' if grid[new_x][new_y] == '@' else '%'
                    new_grid[box_x][box_y] = '@' if grid[box_x][box_y] == '-' else '$'
                    new_players[player_idx] = (new_x, new_y)
                    next_states.append((move, new_grid, new_players))
    
    return next_states

def is_solved(grid):
    goals = sum(row.count('X') for row in grid)
    boxes_on_goals = sum(row.count('$') for row in grid)
    return goals == 0 and boxes_on_goals == 3

def grid_to_string(grid):
    return ''.join(''.join(row) for row in grid)

def solve_sokoban(initial_grid):
    players = []
    for i in range(len(initial_grid)):
        for j in range(len(initial_grid[0])):
            if initial_grid[i][j] in ['*', '%']:
                players.append((i, j))
    
    visited = set()
    queue = deque([([], initial_grid, players)])
    
    while queue:
        moves, grid, current_players = queue.popleft()
        grid_str = grid_to_string(grid)
        
        if grid_str in visited:
            continue
            
        visited.add(grid_str)
        
        if is_solved(grid):
            print(''.join(moves))
            return
            
        for move, new_grid, new_players in get_next_states(grid, current_players):
            queue.append((moves + [move], new_grid, new_players))

# Initialize the puzzle
initial_grid = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', 'X', '@', '*', '+'],
    ['+', '+', '+', '+', '@', '@', '+'],
    ['+', '+', '+', 'X', '-', '-', '+'],
    ['+', '$', '+', '+', '+', 'X', '+'],
    ['+', '+', '+', '$', '+', '+', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solve_sokoban(initial_grid)