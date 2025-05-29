def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def find_positions(grid, chars):
    positions = []
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in chars:
                positions.append((i, j))
    return positions

def evaluate_state(grid):
    boxes = find_positions(grid, ['@', '$'])
    goals = find_positions(grid, ['X', '$', '%'])
    
    # Sum of manhattan distances from each box to nearest goal
    total_distance = 0
    for box in boxes:
        min_dist = min(manhattan_distance(box, goal) for goal in goals)
        total_distance += min_dist
    return total_distance

def dfs_with_limit(grid, depth_limit=30):
    player_pos = None
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    def is_valid_move(x, y):
        return 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != '+'
    
    def try_move(pos, direction, path, depth, visited):
        if depth >= depth_limit:
            return None
        
        px, py = pos
        dx, dy = direction
        new_x, new_y = px + dx, py + dy
        
        if not is_valid_move(new_x, new_y):
            return None
        
        new_grid = [row[:] for row in grid]
        move_char = ''
        
        if direction == (0, 1): move_char = 'R'
        elif direction == (0, -1): move_char = 'L'
        elif direction == (-1, 0): move_char = 'U'
        elif direction == (1, 0): move_char = 'D'
        
        # Moving to empty space or goal
        if grid[new_x][new_y] in ['-', 'X']:
            new_grid[px][py] = '-' if grid[px][py] == '*' else 'X'
            new_grid[new_x][new_y] = '*'
            state_str = str(new_grid)
            if state_str not in visited:
                visited.add(state_str)
                return (new_grid, (new_x, new_y), path + move_char)
        
        # Moving box
        elif grid[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if is_valid_move(box_x, box_y) and grid[box_x][box_y] in ['-', 'X']:
                new_grid[px][py] = '-' if grid[px][py] == '*' else 'X'
                new_grid[new_x][new_y] = '*'
                new_grid[box_x][box_y] = '@' if grid[box_x][box_y] == '-' else '$'
                state_str = str(new_grid)
                if state_str not in visited:
                    visited.add(state_str)
                    return (new_grid, (new_x, new_y), path + move_char)
        
        return None
    
    def is_solved(g):
        return all(c != '@' and c != 'X' for row in g for c in row)
    
    directions = [(0, 1), (0, -1), (-1, 0), (1, 0)]
    stack = [(grid, player_pos, "", 0)]
    visited = {str(grid)}
    
    while stack:
        current_grid, pos, path, depth = stack.pop()
        
        if is_solved(current_grid):
            return path
        
        # Sort moves by heuristic value
        moves = []
        for d in directions:
            result = try_move(pos, d, path, depth, visited)
            if result:
                new_grid, new_pos, new_path = result
                score = evaluate_state(new_grid)
                moves.append((score, (new_grid, new_pos, new_path, depth + 1)))
        
        # Add moves to stack in order of increasing score (better moves last)
        moves.sort(reverse=True)
        for _, move in moves:
            stack.append(move)
    
    return None

# Initial grid
initial_grid = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', 'X', '-', '@', '-', '-', '-', '+', '+'],
    ['+', '+', '+', 'X', '-', '@', '$', '-', 'X', '+'],
    ['+', '+', '-', '-', '@', '-', '$', '-', '-', '+'],
    ['+', '-', '-', '-', '-', 'X', '+', '+', '@', '+'],
    ['+', '-', '-', '-', '-', '@', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', 'X', '+'],
    ['+', '-', '-', '-', '@', '*', '@', 'X', '$', '+'],
    ['+', '+', '+', '-', '-', '-', '-', '-', 'X', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = dfs_with_limit(initial_grid)
print(solution)