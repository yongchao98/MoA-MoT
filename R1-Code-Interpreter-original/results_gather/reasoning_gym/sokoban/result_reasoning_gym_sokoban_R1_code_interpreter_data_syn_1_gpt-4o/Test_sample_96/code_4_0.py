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
    
    # Check if the player can move to the new position
    if grid[x][y] in ('-', 'X'):
        return True
    # Check if the player can push a box
    elif grid[x][y] in ('@', '$'):
        # Calculate the position behind the box
        if move == 'U':
            x -= 1
        elif move == 'D':
            x += 1
        elif move == 'L':
            y -= 1
        elif move == 'R':
            y += 1
        # Check if the position behind the box is empty or a goal
        if grid[x][y] in ('-', 'X'):
            return True
    return False

def move_player_and_box(grid, player_pos, move):
    x, y = player_pos
    if move == 'U':
        x -= 1
    elif move == 'D':
        x += 1
    elif move == 'L':
        y -= 1
    elif move == 'R':
        y += 1
    
    new_grid = [list(row) for row in grid]
    # Check if the player is pushing a box
    if new_grid[x][y] in ('@', '$'):
        # Calculate the position behind the box
        box_x, box_y = x, y
        if move == 'U':
            box_x -= 1
        elif move == 'D':
            box_x += 1
        elif move == 'L':
            box_y -= 1
        elif move == 'R':
            box_y += 1
        # Move the box
        if new_grid[box_x][box_y] == '-':
            new_grid[box_x][box_y] = '@'
        elif new_grid[box_x][box_y] == 'X':
            new_grid[box_x][box_y] = '$'
        # Update the box's previous position
        if new_grid[x][y] == '@':
            new_grid[x][y] = '-'
        elif new_grid[x][y] == '$':
            new_grid[x][y] = 'X'
    
    # Move the player
    if new_grid[x][y] == '-':
        new_grid[x][y] = '%'
    elif new_grid[x][y] == 'X':
        new_grid[x][y] = '%'
    
    # Update the player's previous position
    if new_grid[player_pos[0]][player_pos[1]] == '%':
        new_grid[player_pos[0]][player_pos[1]] = '-'
    elif new_grid[player_pos[0]][player_pos[1]] == '%':
        new_grid[player_pos[0]][player_pos[1]] = 'X'
    
    return new_grid, (x, y)

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
                new_grid, new_player_pos = move_player_and_box(current_grid, current_player_pos, move)
                
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