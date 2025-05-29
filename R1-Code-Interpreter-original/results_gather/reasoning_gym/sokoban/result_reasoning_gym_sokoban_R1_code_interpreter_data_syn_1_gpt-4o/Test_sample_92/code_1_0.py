from collections import deque

def is_valid_position(grid, x, y):
    return grid[x][y] in ['-', 'X']

def move_position(pos, direction):
    x, y = pos
    if direction == 'U':
        return (x - 1, y)
    elif direction == 'D':
        return (x + 1, y)
    elif direction == 'L':
        return (x, y - 1)
    elif direction == 'R':
        return (x, y + 1)

def is_solved(grid):
    for row in grid:
        for cell in row:
            if cell == '@':
                return False
    return True

def sokoban_solver(grid, player_pos):
    directions = ['U', 'D', 'L', 'R']
    queue = deque([(grid, player_pos, "")])
    visited = set()
    visited.add((tuple(map(tuple, grid)), player_pos))

    while queue:
        current_grid, current_player_pos, path = queue.popleft()

        if is_solved(current_grid):
            return path

        for direction in directions:
            new_player_pos = move_position(current_player_pos, direction)
            x, y = new_player_pos

            if is_valid_position(current_grid, x, y):
                new_grid = [list(row) for row in current_grid]
                if new_grid[x][y] == '@':
                    new_box_pos = move_position(new_player_pos, direction)
                    bx, by = new_box_pos
                    if is_valid_position(new_grid, bx, by):
                        new_grid[x][y] = '-'
                        new_grid[bx][by] = '@'
                        new_grid[current_player_pos[0]][current_player_pos[1]] = '-'
                        new_grid[x][y] = '*'
                        if (tuple(map(tuple, new_grid)), new_player_pos) not in visited:
                            visited.add((tuple(map(tuple, new_grid)), new_player_pos))
                            queue.append((new_grid, new_player_pos, path + direction))
                else:
                    new_grid[current_player_pos[0]][current_player_pos[1]] = '-'
                    new_grid[x][y] = '*'
                    if (tuple(map(tuple, new_grid)), new_player_pos) not in visited:
                        visited.add((tuple(map(tuple, new_grid)), new_player_pos))
                        queue.append((new_grid, new_player_pos, path + direction))

    return "No solution"

initial_grid = [
    ['+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', 'X', '@', '-', '$', '+'],
    ['+', '*', '-', '-', '@', '-', '-', '+'],
    ['+', '-', '@', 'X', '-', '-', 'X', '+'],
    ['+', '-', '$', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+']
]

player_start_pos = (2, 1)

solution = sokoban_solver(initial_grid, player_start_pos)
print(solution)