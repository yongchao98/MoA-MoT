from collections import deque
import copy

def is_valid_pos(x, y, width, height):
    return 0 <= x < height and 0 <= y < width

def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def convert_to_matrix(puzzle):
    return [list(row.strip()) for row in puzzle.strip().split('\n')]

def get_state_key(matrix):
    return ''.join([''.join(row) for row in matrix])

def is_goal_state(matrix):
    boxes = 0
    goals = 0
    box_on_goals = 0
    for row in matrix:
        for cell in row:
            if cell in ['@', '$']: boxes += 1
            if cell in ['X', '$', '%', '*']: goals += 1
            if cell in ['$', '*']: box_on_goals += 1
    return boxes == goals == box_on_goals

def solve_sokoban(puzzle):
    # Initialize
    matrix = convert_to_matrix("""
+ + + + + + + + + +
+ - - - X - - - $ +
+ - - @ - @ - - + +
+ - X - - X - @ X +
+ - - - - - - $ @ +
+ - - - - - - $ * +
+ + + + + + + + + +
""")
    
    height = len(matrix)
    width = len(matrix[0])
    
    # Find player position
    player_pos = None
    for i in range(height):
        for j in range(width):
            if matrix[i][j] in ['@', '%', '*']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    # BFS
    queue = deque([(matrix, player_pos, "")])
    visited = {get_state_key(matrix)}
    
    while queue:
        current_matrix, (px, py), path = queue.popleft()
        
        if is_goal_state(current_matrix):
            return path
        
        for direction in ['U', 'D', 'L', 'R']:
            next_px, next_py = get_next_pos(px, py, direction)
            
            if not is_valid_pos(next_px, next_py, width, height):
                continue
                
            if current_matrix[next_px][next_py] == '+':
                continue
                
            new_matrix = copy.deepcopy(current_matrix)
            
            # If next position has a box
            if new_matrix[next_px][next_py] in ['@', '$']:
                box_next_x, box_next_y = get_next_pos(next_px, next_py, direction)
                
                if not is_valid_pos(box_next_x, box_next_y, width, height):
                    continue
                    
                if new_matrix[box_next_x][box_next_y] in ['+', '@', '$']:
                    continue
                
                # Move box
                if new_matrix[box_next_x][box_next_y] in ['X', '%']:
                    new_matrix[box_next_x][box_next_y] = '$'
                else:
                    new_matrix[box_next_x][box_next_y] = '@'
                    
                if new_matrix[next_px][next_py] == '$':
                    new_matrix[next_px][next_py] = 'X'
                else:
                    new_matrix[next_px][next_py] = '-'
                    
            # Move player
            if new_matrix[px][py] in ['*', '%']:
                new_matrix[px][py] = 'X'
            else:
                new_matrix[px][py] = '-'
                
            if new_matrix[next_px][next_py] == 'X':
                new_matrix[next_px][next_py] = '%'
            else:
                new_matrix[next_px][next_py] = '*'
            
            new_state_key = get_state_key(new_matrix)
            if new_state_key not in visited:
                visited.add(new_state_key)
                queue.append((new_matrix, (next_px, next_py), path + direction))
    
    return "No solution found"

# Solve the puzzle
solution = solve_sokoban("")
print(solution)