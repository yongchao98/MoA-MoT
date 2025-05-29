from collections import deque

def is_valid_move(grid, player_pos, move):
    x, y = player_pos
    if move == 'U':
        x -= 1
    elif move == 'D':
        x += 1
    elif move == 'L':
        y -= 1
    elif move == 'R':
        y += 1
    if grid[x][y] in ('-', 'X'):
        return True
    return False

def move_player(grid, player_pos, move):
    x, y = player_pos
    if move == 'U':
        x -= 1
    elif move == 'D':
        x += 1
    elif move == 'L':
        y -= 1
    elif move == 'R':
        y += 1
    return (x, y)

def is_solved(grid):
    for row in grid:
        for cell in row:
            if cell == '@':
                return False
    return True

def solve_sokoban(grid, player_pos):
    moves = ['U', 'D', 'L', 'R']
    queue = deque([(grid, player_pos, "")])
    visited = set()
    visited.add((player_pos, tuple(map(tuple, grid))))
    
    while queue:
        current_grid, current_player_pos, path = queue.popleft()
        
        if is_solved(current_grid):
            return path
        
        for move in moves:
            if is_valid_move(current_grid, current_player_pos, move):
                new_player_pos = move_player(current_grid, current_player_pos, move)
                new_grid = [list(row) for row in current_grid]
                new_grid[current_player_pos[0]][current_player_pos[1]] = '-'
                new_grid[new_player_pos[0]][new_player_pos[1]] = '%'
                
                if (new_player_pos, tuple(map(tuple, new_grid))) not in visited:
                    visited.add((new_player_pos, tuple(map(tuple, new_grid))))
                    queue.append((new_grid, new_player_pos, path + move))
    
    return "No solution"

# Initial grid setup
grid = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '$', '%', '@', '-', '+'],
    ['+', '+', '+', '@', '$', '-', '+'],
    ['+', '$', '+', '-', '-', '-', '+'],
    ['+', '+', '-', '-', '-', '-', '+'],
    ['+', '+', '-', '-', '-', '-', '+'],
    ['+', '+', '-', '-', '@', 'X', '+'],
    ['+', 'X', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

player_pos = (1, 3)

solution = solve_sokoban(grid, player_pos)
print(solution)